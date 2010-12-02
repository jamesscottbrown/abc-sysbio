# cuda utilities

import os, time, string

import pycuda.driver as cuda
import pycuda.tools as tools
from pycuda.compiler import SourceModule
import pycuda.autoinit
import numpy as np
from struct import unpack

_lsoda_source_ = """

extern "C"{

__device__ myFex myfex;
__device__ myJex myjex;

__global__ void init_common(){
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  cuLsodaCommonBlockInit( &(common[tid]) );
}

__global__ void cuLsoda(int *neq, double *y, double *t, double *tout, int *itol, 
			double *rtol, double *atol, int *itask, int *istate, int *iopt, 
                        double *rwork, int *lrw, int *iwork, int *liw, int *jt)
{
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  
  dlsoda_(myfex, neq+tid, y+tid*NXMAX, t+tid, tout+tid, itol+tid, rtol+tid, atol+tid, itask+tid, 
	  istate+tid, iopt+tid, rwork+tid*RSIZE, lrw+tid, iwork+tid*ISIZE, liw+tid, myjex, jt+tid, &(common[tid]) );
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

def compile_lsoda( application_code, cuLsoda, blockSize, gridSize, pmax, xmax, options, write=False ):

    fa = open(application_code,'r')
    fc = open(cuLsoda,'r')

    isize = 20 + xmax
    rsize = 22 + xmax * max(16, xmax + 9)

    _isize_ = "#define ISIZE " + repr( 20 + xmax ) + "\n"
    _rsize_ = "#define RSIZE " + repr( 22 + xmax * max(16, xmax + 9) ) + "\n"

    _xmax_ = "#define NXMAX " + repr(xmax) + "\n"
    _textures_ = "texture<int, 1, cudaReadModeElementType> model_tex;\n" + "texture<float, 2, cudaReadModeElementType> param_tex;\n"
    _common_block_ = "__device__ struct cuLsodaCommonBlock common[" + repr(blockSize*gridSize) + "];\n"

    _code_ =  _isize_ + _rsize_ + _xmax_ + _textures_ + fa.read() + fc.read() + _common_block_ + _lsoda_source_

    if write==True:
        fo = open('source_o.cu','w')
        print >>fo, _code_

    compiled = pycuda.compiler.SourceModule( _code_, nvcc="nvcc", options=options, no_extern_c=True )

    model_tex = compiled.get_texref("model_tex")
    param_tex = compiled.get_texref("param_tex")
    init_common_Kernel = compiled.get_function("init_common")
    init_common_Kernel( block=(blockSize,1,1), grid=(gridSize,1) )

    lsoda_Kernel = compiled.get_function("cuLsoda")

    print "#### max threads / registers / blocks :", tools.DeviceData().max_threads, tools.DeviceData().registers, tools.DeviceData().thread_blocks_per_mp
    print "#### kernel mem local / shared / registers : ", lsoda_Kernel.local_size_bytes, lsoda_Kernel.shared_size_bytes, lsoda_Kernel.num_regs
    occ = tools.OccupancyRecord( tools.DeviceData(), threads=blockSize, shared_mem=lsoda_Kernel.shared_size_bytes, registers=lsoda_Kernel.num_regs )
    print "#### threadblocks per mp / limit / occupancy :", occ.tb_per_mp, occ.limited_by, occ.occupancy

    return [ lsoda_Kernel, model_tex, param_tex ]
        
def run_lsoda(ptrs, blockSize, gridSize, nm, pmax, xmax, neqn, model, param, times, init, dt, in_atol=1e-12, in_rtol=1e-6 ):
    start_time = time.time()

    lsoda_Kernel = ptrs[0]
    model_tex = ptrs[1]
    param_tex = ptrs[2]

    nt = blockSize * gridSize

    # output array
    ntimes = len(times)
    ret_xt = np.zeros( [nt, ntimes, xmax] )

    # calculate sizes of work spaces
    isize = 20 + xmax
    rsize = 22 + xmax * max(16, xmax + 9)
        
    # local variables
    t      = np.zeros( [nt], dtype=np.float64)
    jt     = np.zeros( [nt], dtype=np.int32)
    neq    = np.zeros( [nt], dtype=np.int32)
    itol   = np.zeros( [nt], dtype=np.int32)
    iopt   = np.zeros( [nt], dtype=np.int32)
    rtol   = np.zeros( [nt], dtype=np.float64)
    iout   = np.zeros( [nt], dtype=np.int32)
    tout   = np.zeros( [nt], dtype=np.float64)
    itask  = np.zeros( [nt], dtype=np.int32)
    istate = np.zeros( [nt], dtype=np.int32)
    atol   = np.zeros( [nt], dtype=np.float64)

    liw    = np.zeros( [nt], dtype=np.int32)
    lrw    = np.zeros( [nt], dtype=np.int32)
    iwork  = np.zeros( [isize*nt], dtype=np.int32)
    rwork  = np.zeros( [rsize*nt], dtype=np.float64)
    y      = np.zeros( [xmax*nt], dtype=np.float64)

    # initialise variables
    for i in range(nt):
        neq[i] = neqn[model[i]]
        t[i] = times[0]
        itol[i] = 1
        itask[i] = 1
        istate[i] = 1
        iopt[i] = 0
        jt[i] = 2
        atol[i] = in_atol
        rtol[i] = in_rtol

        liw[i] = isize
        lrw[i] = rsize

        # initial conditions
        for im in range(nm):
            if model[i] == im:
                for k in range(len(init[im])):
                    y[i*xmax + k] = init[im][k]
                    ret_xt[i, 0, k] = init[im][k]

    # allocate on device
    d_t      = cuda.mem_alloc(t.size      * t.dtype.itemsize)
    d_jt     = cuda.mem_alloc(jt.size     * jt.dtype.itemsize)
    d_neq    = cuda.mem_alloc(neq.size    * neq.dtype.itemsize)
    d_liw    = cuda.mem_alloc(liw.size    * liw.dtype.itemsize)
    d_lrw    = cuda.mem_alloc(lrw.size    * lrw.dtype.itemsize)
    d_itol   = cuda.mem_alloc(itol.size   * itol.dtype.itemsize)
    d_iopt   = cuda.mem_alloc(iopt.size   * iopt.dtype.itemsize)
    d_rtol   = cuda.mem_alloc(rtol.size   * rtol.dtype.itemsize)
    d_iout   = cuda.mem_alloc(iout.size   * iout.dtype.itemsize)
    d_tout   = cuda.mem_alloc(tout.size   * tout.dtype.itemsize)
    d_itask  = cuda.mem_alloc(itask.size  * itask.dtype.itemsize)
    d_istate = cuda.mem_alloc(istate.size * istate.dtype.itemsize)
    d_y      = cuda.mem_alloc(y.size      * y.dtype.itemsize)
    d_atol   = cuda.mem_alloc(atol.size   * atol.dtype.itemsize)
    d_iwork  = cuda.mem_alloc(iwork.size  * iwork.dtype.itemsize)
    d_rwork  = cuda.mem_alloc(rwork.size  * rwork.dtype.itemsize)

    # copy to device
    cuda.memcpy_htod(d_t, t)
    cuda.memcpy_htod(d_jt, jt)
    cuda.memcpy_htod(d_neq, neq)
    cuda.memcpy_htod(d_liw, liw)
    cuda.memcpy_htod(d_lrw, lrw)
    cuda.memcpy_htod(d_itol, itol)
    cuda.memcpy_htod(d_iopt, iopt)
    cuda.memcpy_htod(d_rtol, rtol)
    cuda.memcpy_htod(d_iout, iout)
    cuda.memcpy_htod(d_tout, tout)
    cuda.memcpy_htod(d_itask, itask)
    cuda.memcpy_htod(d_istate, istate)
    cuda.memcpy_htod(d_y, y)
    cuda.memcpy_htod(d_atol, atol)
    cuda.memcpy_htod(d_iwork, iwork)
    cuda.memcpy_htod(d_rwork, rwork)

    # model texture
    d_m = cuda.mem_alloc(model.size*model.dtype.itemsize)
    cuda.memcpy_htod(d_m, model.astype(np.int32) )
    model_tex.set_address(d_m, model.size*model.dtype.itemsize, allow_offset=False)

    # parameter texture
    ary = create_2D_array( param )
    copy2D_host_to_array(ary, param, pmax*4, nt )
    param_tex.set_array(ary)

    #print "Running lsoda"

    if dt < 0:
        start_time = time.time()
        for i in range(1,len(times)):

            for j in range(nt):
                tout[j] = times[i]; 
            cuda.memcpy_htod( d_tout, tout ) 

            lsoda_Kernel( d_neq, d_y, d_t, d_tout, d_itol, d_rtol, d_atol, d_itask, d_istate,
                          d_iopt, d_rwork, d_lrw, d_iwork, d_liw, d_jt, block=(blockSize,1,1), grid=(gridSize,1) );

            cuda.memcpy_dtoh(t, d_t)
            cuda.memcpy_dtoh(y, d_y)
            cuda.memcpy_dtoh(istate, d_istate)

            for j in range(nt):
                for k in range(xmax):
                    ret_xt[j, i, k] = y[j*xmax + k]

        # end of loop over time points

    else:
        tt = times[0]
        
        start_time = time.time()
        for i in range(1,len(times)):
            
            while 1:
                
                next_time = min(tt+dt, times[i])
                
                for j in range(nt):
                    tout[j] = next_time; 
                cuda.memcpy_htod( d_tout, tout ) 

                lsoda_Kernel( d_neq, d_y, d_t, d_tout, d_itol, d_rtol, d_atol, d_itask, d_istate,
                              d_iopt, d_rwork, d_lrw, d_iwork, d_liw, d_jt, block=(blockSize,1,1), grid=(gridSize,1) );

                cuda.memcpy_dtoh(t, d_t)
                cuda.memcpy_dtoh(y, d_y)
                cuda.memcpy_dtoh(istate, d_istate)

                if np.abs(next_time - times[i]) < 1e-5:
                    tt = next_time
                    break

                tt = next_time
                

            for j in range(nt):
                for k in range(xmax):
                    ret_xt[j, i, k] = y[j*xmax + k]

        # end of loop over time points

    #print "### ###"
    #print istate
    return ret_xt


    



