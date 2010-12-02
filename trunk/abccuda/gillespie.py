# cuda utilities

import os, time, string

import pycuda.driver as cuda
import pycuda.tools as tools
from pycuda.compiler import SourceModule
import pycuda.autoinit
import numpy as np
from struct import unpack

# total device global memory usage
# MT RNG
# MT  : 32768 * sizeof(mt_struct_stripped) = 32768 x 16 = 524288 bytes
# MTS : 32768 * sizeof(MersenneTwisterState) = 32768 x 96 = 3145728 bytes
# values
# t : nt x 4
# x : nt x pitch (usually 64)
# tot = 3670016 + nt x 68 
# nt = 2048 mem = 3,801,088 = 3.8MB

_rng_test_source_ ="""

__global__ void TestMersenneTwisters(float *z, int *r) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  int nrand = r[tid];

  for(int i=0; i<nrand; ++i){
    z[ i + tid*nrand ] = MersenneTwisterGenerate(&(MTS[tid]), tid)/ 4294967295.0f;
  }
}

"""

_gillespie_source_ ="""

__constant__ float tmax;
__constant__ int vxp;
__constant__ int vxtp;

__device__ int sample(int nh, float* h, float u){
  
  int i = 0;
  for(i=0; i<nh; ++i){
    if( u < h[i] ) break;
    u = u - h[i];
  }
  return i;
}

__global__ void Gillespie_one_step(float *vt, int* vx, int* vxt){

  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  int i, x[NXMAX];

  float h0, h[NRMAX];
  for(i=0; i<NRMAX; ++i){
    h[i] = 0;
  }
 
  for(i=0; i<NXMAX; ++i){
    x[i] = ((int*)( (char*) vx + tid * vxp))[i];
  }
   
  float t  = vt[ tid ];
   
  while( t < tmax ){

    // rules
    rules( x, t );

    // calculate hazards
    hazards( x, h, tid );

    h0 = 0;
    for(i=0; i<NRMAX; ++i) h0 += h[i];

    if(h0 <= 0){
      
      // copy the current local values over to xt
      for(i=0; i<NXMAX; ++i){
        ((int*)( (char*) vxt + tid * vxtp))[i] = x[i];
      }
      break;
      
    }else{
      float u1 = MersenneTwisterGenerate(&(MTS[tid]), tid)/ 4294967295.0f;
      float u2 = MersenneTwisterGenerate(&(MTS[tid]), tid)/ 4294967295.0f;

      // increment the time
      t += -log(u1)/h0;
    
      if( t >= tmax ){

        // copy the current local values over to xt
        for(i=0; i<NXMAX; ++i){
          ((int*)( (char*) vxt + tid * vxtp))[i] = x[i];
        }  
      }

      // sample reaction
      i = sample(NRMAX, h, h0*u2);

      // update stochiometry
      stoichiometry( x ,i, tid );
      }
    }

    vt[tid] = t;
    for(i=0; i<NXMAX; ++i){
      ((int*)( (char*) vx + tid * vxp))[i] = x[i];
    }
}

"""


def copy2D_host_to_device(dev, host, src_pitch, dst_pitch, width, height ):
    copy = cuda.Memcpy2D()
    copy.set_src_host(host)
    copy.set_dst_device(dev)
    copy.src_pitch = src_pitch
    copy.dst_pitch = dst_pitch
    copy.width_in_bytes = width
    copy.height = height
    copy(aligned=True)

def copy2D_device_to_host(host, dev, src_pitch, dst_pitch, width, height ):
    copy = cuda.Memcpy2D()
    copy.set_src_device(dev)
    copy.set_dst_host(host)
    copy.src_pitch = src_pitch
    copy.dst_pitch = dst_pitch
    copy.width_in_bytes = width
    copy.height = height
    copy(aligned=True)

# Create a 2D GPU array (for assignment
# to texture) from a numpy 2D array
def create_2D_array( mat ):
    descr = cuda.ArrayDescriptor()
    descr.width = mat.shape[1]
    descr.height = mat.shape[0]
    descr.format = cuda.dtype_to_array_format( mat.dtype )
    descr.num_channels = 1
    descr.flags = 0
    ary = cuda.Array(descr)
    return ary

# Copy 2D host numpy array to 2D
# GPU array object
def copy2D_host_to_array(arr, host, width, height ):
    copy = cuda.Memcpy2D()
    copy.set_src_host(host)
    copy.set_dst_array(arr)
    copy.height = height
    copy.width_in_bytes = copy.src_pitch = width
    copy.height = height
    copy(aligned=True)

def compile_mt_code( code, mt_cu, mt_data, blockSize, gridSize, options ):

    # get the location of the MersenneTwister.cu file
    #s = abcsysbio.__file__
    #mt_cu = s.replace('__init__.pyc','MersenneTwister.cu')
    #print '#### using MersenneTwister.cu at', mt_cu

    f = open(mt_cu,'r')
    _code_ = f.read() + code

    opts = pycuda.driver.jit_option()
    compiled = pycuda.compiler.SourceModule( _code_, nvcc="nvcc", options=options )
    initialise_twisters( mt_data, compiled, blockSize, gridSize )

    return [ compiled ]

def initialise_twisters( mt_data, mod, blockSize, gridSize ):

    # get the location of the MersenneTwister.dat file
    #s = abcsysbio.__file__
    #mt_data = s.replace('__init__.pyc','MersenneTwister.dat')
    #print '#### using MersenneTwister.dat at', mt_data
    
    pMT = int(mod.get_global("MT")[0])
    
    MT_RNG_COUNT = 32768

    f = open(mt_data, 'rb')
    s = f.read(16*MT_RNG_COUNT)
    # the file should contain 32768 x 4 integers
    tup = unpack('131072i',s)

    # Copy the offline MT parameters over to GPU
    cuda.memcpy_htod( pMT, np.array(tup).astype(np.uint32) )

    # Seed the MT
    for i in range(MT_RNG_COUNT):
        place = i*16 + 12
        cuda.memcpy_htod( pMT + place, np.uint32( np.random.randint(100000000)) )

    InitialiseAllMersenneTwisters = mod.get_function("InitialiseAllMersenneTwisters")
    InitialiseAllMersenneTwisters( block=(blockSize,1,1), grid=(gridSize,1) )

def run_mt_test( mods, blockSize, gridSize, nrandom ):
    nr = np.array([nrandom for i in range(blockSize*gridSize)], dtype=np.int32)
    d_nr = cuda.mem_alloc(nr.size*nr.dtype.itemsize)
    cuda.memcpy_htod(d_nr, nr)

    r = np.zeros( blockSize*gridSize*nrandom, np.float32)
    d_r = cuda.mem_alloc(r.size*r.dtype.itemsize)
    cuda.memcpy_htod(d_r, r)

    TestMersenneTwisters = mods[0].get_function("TestMersenneTwisters")
    TestMersenneTwisters( d_r, d_nr, block=(blockSize,1,1), grid=(gridSize,1) )

    cuda.memcpy_dtoh(r, d_r)

    return r

def compile_gillespie( application_code, mt_cu, mt_data, blockSize, gridSize, pmax, xmax, options, write=False ):

    f = open(application_code,'r')
    _source_ = '#define NPMAX ' + repr(pmax) + '\n' + '#define NXMAX ' + repr(xmax) + '\n' + 'texture<int, 1, cudaReadModeElementType> model_tex;\n' + 'texture<float, 2, cudaReadModeElementType> param_tex;\n' + f.read() +_gillespie_source_

    m = compile_mt_code( _source_, mt_cu, mt_data, blockSize, gridSize, options)

    if write==True:
        fo = open('source_g.cu','w')
        print >>fo, _source_

    ptmax = int(m[0].get_global("tmax")[0])
    pvxp = int(m[0].get_global("vxp")[0])
    pvxtp = int(m[0].get_global("vxtp")[0])
    model_tex = m[0].get_texref("model_tex")
    param_tex = m[0].get_texref("param_tex")
    Gillespie_Kernel = m[0].get_function("Gillespie_one_step")

    print "#### max threads / registers / blocks :", tools.DeviceData().max_threads, tools.DeviceData().registers, tools.DeviceData().thread_blocks_per_mp
    print "#### kernel mem local / shared / registers : ", Gillespie_Kernel.local_size_bytes, Gillespie_Kernel.shared_size_bytes, Gillespie_Kernel.num_regs
    occ = tools.OccupancyRecord( tools.DeviceData(), threads=blockSize, shared_mem=Gillespie_Kernel.shared_size_bytes, registers=Gillespie_Kernel.num_regs )
    print "#### threadblocks per mp / limit / occupancy :", occ.tb_per_mp, occ.limited_by, occ.occupancy

    return [ Gillespie_Kernel, ptmax, pvxp, pvxtp, model_tex, param_tex ]
        
def run_gillespie(ptrs, blockSize, gridSize, nm, pmax, xmax, model, param, times, init ):
    start_time = time.time()

    Gillespie_Kernel = ptrs[0]
    ptmax = ptrs[1]
    pvxp = ptrs[2]
    pvxtp = ptrs[3]
    model_tex = ptrs[4]
    param_tex = ptrs[5]

    nt = blockSize * gridSize

    # output array
    ntimes = len(times)
    ret_xt = np.zeros( [nt, ntimes, xmax] )

    # times array
    t = np.zeros( [nt] ).astype(np.float32)
    d_t = cuda.mem_alloc(t.size*t.dtype.itemsize)
    #print '#### times array:', time.time() - start_time

    # 2D species arrays
    x = np.zeros( [nt,xmax], dtype=np.int32)
    d_x , p_x = cuda.mem_alloc_pitch(width=xmax*4, height=nt, access_size=4 )
    xt = np.zeros( [nt,xmax], dtype=np.int32)
    d_xt , p_xt = cuda.mem_alloc_pitch(width=xmax*4, height=nt, access_size=4 )
    #print '#### 2D arrays:', time.time() - start_time

    cuda.memcpy_htod(pvxp, np.array([p_x], dtype=np.int32) )
    cuda.memcpy_htod(pvxtp, np.array([p_xt], dtype=np.int32) )

    # initialize all arrays
    for ip in range(nt):
        t[ip] = times[0]

        for im in range(nm):
            if model[ip] == im:
                for k in range(len(init[im])):
                    x[ip,k] = init[im][k]
                    xt[ip,k] = init[im][k]
                    ret_xt[ip,0,k] = init[im][k]

    # copy to device
    cuda.memcpy_htod(d_t, t)
    copy2D_host_to_device(d_x, x, xmax*4, p_x, xmax*4, nt )
    copy2D_host_to_device(d_xt, xt, xmax*4, p_xt, xmax*4, nt )

    # model texture
    d_m = cuda.mem_alloc(model.size*model.dtype.itemsize)
    cuda.memcpy_htod(d_m, model.astype(np.int32) )
    model_tex.set_address(d_m, model.size*model.dtype.itemsize, allow_offset=False)

    # parameter texture
    ary = create_2D_array( param )
    copy2D_host_to_array(ary, param, pmax*4, nt )
    param_tex.set_array(ary)

    #start_time = time.time()
    for i in range(1,len(times)):
        
        cuda.memcpy_htod(ptmax, np.array([times[i]], dtype=np.float32 ))
        
        Gillespie_Kernel( d_t, d_x, d_xt, block=(blockSize,1,1), grid=(gridSize,1) );
    
        copy2D_device_to_host(xt, d_xt, p_xt, xmax*4, xmax*4, nt )
        
        for j in range(nt):
            for k in range(xmax):
                ret_xt[j,i,k] = xt[j,k]

    return ret_xt




    



