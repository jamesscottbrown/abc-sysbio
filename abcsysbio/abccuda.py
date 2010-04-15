# cuda utilities

import os, string, abcsysbio

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

_mersenne_twister_source_ = """
    // This mersenne twister code is derived from the code at 
    // http://www.jcornwall.me.uk/2009/04/mersenne-twisters-in-cuda/
    // which in turn is derived from the NVIDIA CUDA SDK example 'MersenneTwister'.

    #define MT_RNG_COUNT 32768
    #define MT_MM     9
    #define MT_NN     19
    #define MT_WMASK  0xFFFFFFFFU
    #define MT_UMASK  0xFFFFFFFEU
    #define MT_LMASK  0x1U
    #define MT_RNG_COUNT 32768
    #define MT_SHIFT0 12
    #define MT_SHIFTB 7
    #define MT_SHIFTC 15
    #define MT_SHIFT1 18


    // size = 16 bytes
    struct __align__(16) mt_struct_stripped {
    unsigned int matrix_a;
    unsigned int mask_b;
    unsigned int mask_c;
    unsigned int seed;
    };

    // size = 84 bytes but aligned = 96 bytes
    struct __align__(16) MersenneTwisterState {
    unsigned int mt[MT_NN];
    int iState;
    unsigned int mti1;
    };

    __device__ static mt_struct_stripped MT[MT_RNG_COUNT];
    __device__ static MersenneTwisterState MTS[MT_RNG_COUNT];

    __device__ unsigned int MersenneTwisterGenerate(MersenneTwisterState &state, const unsigned int& threadID) {
    int iState1 = state.iState + 1;
    int iStateM = state.iState + MT_MM;

    if(iState1 >= MT_NN) iState1 -= MT_NN;
    if(iStateM >= MT_NN) iStateM -= MT_NN;

    unsigned int mti = state.mti1;
    state.mti1 = state.mt[iState1];
    unsigned int mtiM = state.mt[iStateM];

    unsigned int x = (mti & MT_UMASK) | (state.mti1 & MT_LMASK);
    x = mtiM ^ (x >> 1) ^ ((x & 1) ? MT[threadID].matrix_a : 0);
    state.mt[state.iState] = x;
    state.iState = iState1;
    
    // Tempering transformation.
    x ^= (x >> MT_SHIFT0);
    x ^= (x << MT_SHIFTB) & MT[threadID].mask_b;
    x ^= (x << MT_SHIFTC) & MT[threadID].mask_c;
    x ^= (x >> MT_SHIFT1);
  
    return x;
    }

    __device__ float runif(float a, float b, MersenneTwisterState &state, const unsigned int& threadID) {
    //Convert to (0, 1] float
    float u = float(MersenneTwisterGenerate(state, threadID)) / 4294967295.0f;
    //Convert to (a,b] float
    return a + (b-a)*u;
    }

    __device__ float rexp(float lambda, MersenneTwisterState &state, const unsigned int& threadID) {
    return -log( runif(0, 1, state, threadID) )/lambda;
    }

    __device__ void MersenneTwisterInitialise(MersenneTwisterState &state, const unsigned int& threadID) {
    state.mt[0] = MT[threadID].seed;
    for(int i = 1; i < MT_NN; ++ i) {
    state.mt[i] = (1812433253U * (state.mt[i - 1] ^ (state.mt[i - 1] >> 30)) + i) & MT_WMASK;
    }
    state.iState = 0;
    state.mti1 = state.mt[0];
    }

    __global__ void InitialiseAllMersenneTwisters() {
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    MersenneTwisterInitialise(MTS[tid],tid);
    //MersenneTwisterBurnIn(MTS[tid], tid);
    }

    __global__ void TestMersenneTwistersUniform(float *z, int *r) {
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int n = r[tid];
    for(int i=0; i<n; ++i){
    z[ i + tid*n ] = runif(0, 1, MTS[tid], tid);
    }
    }

    __global__ void TestMersenneTwistersExponential(float *z, int *r) {
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int n = r[tid];
    for(int i=0; i<n; ++i){
    z[ i + tid*n ] = rexp(1, MTS[tid], tid);
    }
    }
    """

_textures_="""

    texture<int, 1, cudaReadModeElementType> model_tex;
    texture<float, 2, cudaReadModeElementType> param_tex;

    __global__ void TestTextureAccess(int *m, float *p, int* pitch) {
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;

    m[tid] = tex1Dfetch(model_tex, tid);

    int vxp = pitch[tid];
    for(int i=0; i<5; ++i){
       ((float*)( (char*) p + tid * vxp))[i] = tex2D(param_tex, i, tid);
    }

    
    }
"""

_sample_="""
    __device__ int sample(int& nh, float& h0, float* h, MersenneTwisterState &state, const unsigned int& threadID){
         
         float u = runif(0, h0, state, threadID);
         int i = 0;

         for(i=0; i<nh; ++i){
             if( u < h[i] ) break;
             u = u - h[i];
         }
         return i;
    }

    __global__ void TestSample(int *z, int *r){
         const int tid = blockDim.x * blockIdx.x + threadIdx.x;

         int nh = 3;
         float h[] = {0.0,0.0,0.0};

         h[0] = 0.1;
         h[1] = 0.4;
         h[2] = 0.5;
         float h0 = h[0] + h[1] + h[2];
       
         int n = r[tid];
         for(int i=0; i<n; ++i){
             z[ i + tid*n ] = sample(nh,h0,h,MTS[tid],tid);
         }
    }

"""

_gillespie_source_="""
    __constant__ float tmax;
    __constant__ int vxp;

    __global__ void Gillespie_one_step(float* vt, float* vx) {
    const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

    int i;
    int model = tex1Dfetch(model_tex, tid);

    float x[ NXMAX ];
    for(i=0; i<nx[model]; ++i){
       x[i] = ((float*)( (char*) vx + tid * vxp))[i];
    }
    
    float t  = vt[ tid ];
    float h0;
    
    while( t < tmax ){
       hazards(h0, x, model, MTS[tid], tid);
       if(h0 <= 0.0) break;
       t = t + rexp( h0, MTS[tid], tid );
    }

    vt[tid] = t;
    for(i=0; i<nx[model]; ++i){
       ((float*)( (char*) vx + tid * vxp))[i] = x[i];
    }
    }
    """

_test_hazards_source_="""
    __constant__ int vxp;

    __global__ void test_hazards(int* vh, float* vx){
    const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

    int i;
    int model = tex1Dfetch(model_tex, tid);

    float x[ NXMAX ];
    for(i=0; i<nx[model]; ++i){
       x[i] = ((float*)( (char*) vx + tid * vxp))[i];
    }

    float h0;
    int r = hazards(h0, x, model, MTS[tid], tid);
    
    for(i=0; i<nx[model]; ++i){
       ((float*)( (char*) vx + tid * vxp))[i] = x[i];
    }

    vh[tid] = r;
    }
    """

_gillespie2_source_="""
    __constant__ int vxp;

    __global__ void Gillespie_one_step(float* vt, float* vx, float* vtmax) {
    const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

    int i;
    int model = tex1Dfetch(model_tex, tid);

    float tmax = vtmax[ tid ];

    float x[ NXMAX ];
    for(i=0; i<nx[model]; ++i){
       x[i] = ((float*)( (char*) vx + tid * vxp))[i];
    }
    
    float t  = vt[ tid ];
    float h0;
    
    while( t < tmax ){
       hazards(h0, x, model, MTS[tid], tid);
       if(h0 <= 0.0) break;
       t = t + rexp( h0, MTS[tid], tid );
    }

    vt[tid] = t;
    for(i=0; i<nx[model]; ++i){
       ((float*)( (char*) vx + tid * vxp))[i] = x[i];
    }
    }
    """


def initialise_twisters( mt_data, mod, nthread, nblock ):

    # get the location of the MersenneTwister.dat file
    s = abcsysbio.__file__
    mt_data = s.replace('__init__.pyc','MersenneTwister.dat')
    print '#### using MersenneTwister data at', mt_data
    
    print "\n#### initialise_twisters"
    print "#### obtaining MT  pointer of size", mod.get_global("MT")[1], "bytes" 
    print "#### obtaining MTS pointer of size", mod.get_global("MTS")[1], "bytes" 

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
    InitialiseAllMersenneTwisters( block=(nthread,1,1), grid=(nblock,1) )

def compile_mersenne_twister_test( mt_data, nr, nthread, nblock ):

    print "\n#### GPU info"
    print "#### max threads / registers / blocks :", tools.DeviceData().max_threads, tools.DeviceData().registers, tools.DeviceData().thread_blocks_per_mp

    _code_ = _mersenne_twister_source_
    compiled = pycuda.compiler.SourceModule( _code_ )
    initialise_twisters( mt_data, compiled, nthread, nblock )

    twister_test_u = compiled.get_function("TestMersenneTwistersUniform")
    twister_test_e = compiled.get_function("TestMersenneTwistersExponential")

    return [twister_test_u, twister_test_e]

def run_mersenne_twister_test( ptr, nr, nthread, nblock ):
    nrandom = nr*nthread*nblock

    a = np.zeros(nrandom).astype(np.float32)
    d_a = cuda.mem_alloc(a.size*a.dtype.itemsize)
    cuda.memcpy_htod(d_a, a)

    r = (np.ones(nthread*nblock)*nr).astype(np.int32)
    d_r = cuda.mem_alloc(r.size*r.dtype.itemsize)
    cuda.memcpy_htod(d_r, r)
    
    ptr( d_a, d_r, block=(nthread,1,1), grid=(nblock,1) )
    cuda.memcpy_dtoh(a, d_a)

    return a

def compile_sample_test( mt_data, nr, nthread, nblock ):

    print "\n#### GPU info"
    print "#### max threads / registers / blocks :", tools.DeviceData().max_threads, tools.DeviceData().registers, tools.DeviceData().thread_blocks_per_mp

    _code_ = _mersenne_twister_source_ + _sample_
    compiled = pycuda.compiler.SourceModule( _code_ )
    initialise_twisters( mt_data, compiled, nthread, nblock )

    sample_test = compiled.get_function("TestSample")

    return sample_test

def run_sample_test( ptr, nr, nthread, nblock ):
    nrandom = nr*nthread*nblock

    a = np.zeros(nrandom).astype(np.int32)
    d_a = cuda.mem_alloc(a.size*a.dtype.itemsize)
    cuda.memcpy_htod(d_a, a)

    r = (np.ones(nthread*nblock)*nr).astype(np.int32)
    d_r = cuda.mem_alloc(r.size*r.dtype.itemsize)
    cuda.memcpy_htod(d_r, r)
    
    ptr( d_a, d_r, block=(nthread,1,1), grid=(nblock,1) )
    cuda.memcpy_dtoh(a, d_a)

    return a

def compile_texture_test():
    _code_ = _textures_
    compiled = pycuda.compiler.SourceModule( _code_ )
    tex_test = compiled.get_function("TestTextureAccess")
    model_tex = compiled.get_texref("model_tex")
    param_tex = compiled.get_texref("param_tex")

    return [ tex_test, model_tex, param_tex ]

def run_texture_test(ptrs, nthread, nblock, m_in, p_in):

    fun, model_tex, param_tex = ptrs[0], ptrs[1], ptrs[2]

    nt = nblock*nthread
    pmax = 5

    # input model 1D array
    d_m = cuda.mem_alloc(m_in.size*m_in.dtype.itemsize)
    cuda.memcpy_htod(d_m, m_in.astype(np.int32) )
    model_tex.set_address(d_m, m_in.size*m_in.dtype.itemsize, allow_offset=False)

    # input parameters 2D array 
    descr = cuda.ArrayDescriptor()
    descr.width = pmax
    descr.height = nt
    descr.format = cuda.dtype_to_array_format(p_in.dtype)
    descr.num_channels = 1
    descr.flags = 0

    ary = cuda.Array(descr)

    copy = cuda.Memcpy2D()
    copy.set_src_host(p_in)
    copy.set_dst_array(ary)
    copy.height = nt
    copy.width_in_bytes = copy.src_pitch = pmax*4
    copy(aligned=True)

    param_tex.set_array(ary)

    # output model 1D array
    a = np.zeros(nt).astype(np.int32)
    d_a = cuda.mem_alloc(a.size*a.dtype.itemsize)
    cuda.memcpy_htod(d_a, a)

    # output model 2D array
    b = np.zeros( [nt,pmax], order="C", dtype=np.float32)
    vxw = pmax*4
    alloc_tup = cuda.mem_alloc_pitch(width=vxw, height=nt, access_size=4 )
    vxp = alloc_tup[1]
    d_b = alloc_tup[0]

    copy = cuda.Memcpy2D()
    copy.set_src_host(b)
    copy.set_dst_device(d_b)
    copy.src_pitch = vxw
    copy.dst_pitch = vxp
    copy.width_in_bytes = vxw
    copy.height = nt
    copy(aligned=True)

    # output model 2D array pitch
    pitch = ( np.ones( [nt] )*vxp ).astype(np.int32)
    d_pitch = cuda.mem_alloc(pitch.size*pitch.dtype.itemsize)
    cuda.memcpy_htod(d_pitch, pitch)

    # run the kernel
    fun(d_a, d_b, d_pitch, block=(nthread,1,1), grid=(nblock,1) )

    # copy the relevant results back
    cuda.memcpy_dtoh(a, d_a)

    copy = cuda.Memcpy2D()
    copy.set_src_device(d_b)
    copy.set_dst_host(b)
    copy.src_pitch = vxp
    copy.dst_pitch = vxw
    copy.width_in_bytes = vxw
    copy.height = nt
    copy(aligned=True)
    
    return [a, b]
    
def compile_cuda_gillespie( application_code, mt_data, nthread, nblock ):

    print "\n#### GPU info"
    print "#### max threads / registers / blocks :", tools.DeviceData().max_threads, tools.DeviceData().registers, tools.DeviceData().thread_blocks_per_mp

    f = open(application_code,'r')

    _code_ = _mersenne_twister_source_ + _textures_ + _sample_ + f.read() + _gillespie_source_
    
    compiled = pycuda.compiler.SourceModule( _code_ )
    ptmax = int(compiled.get_global("tmax")[0])
    pvxp = int(compiled.get_global("vxp")[0])

    #pnx = int(compiled.get_global("nx")[0])
    #cuda.memcpy_htod(pnx, np.array( [1] ).astype(np.float32) )
    
    print "\n#### compile_cuda_gillespie"
    print "#### obtaining tmax pointer of size", compiled.get_global("tmax")[1], "bytes, value :", compiled.get_global("tmax")[0] 
    print "#### obtaining vxp pointer of size", compiled.get_global("vxp")[1], "bytes, value :", compiled.get_global("vxp")[0]  

    initialise_twisters( mt_data, compiled, nthread, nblock )

    Gillespie_one_step = compiled.get_function("Gillespie_one_step")
    model_tex = compiled.get_texref("model_tex")
    param_tex = compiled.get_texref("param_tex")

    print "\n#### kernel info" 
    print "#### kernel mem local / shared / registers : ", Gillespie_one_step.local_size_bytes, Gillespie_one_step.shared_size_bytes, Gillespie_one_step.num_regs

    print "\n#### occupancy info"
    occ = tools.OccupancyRecord( tools.DeviceData(), threads=nthread, shared_mem=Gillespie_one_step.shared_size_bytes, registers=Gillespie_one_step.num_regs )
    print "#### threadblocks per mp / limit / occupancy :", occ.tb_per_mp, occ.limited_by, occ.occupancy 
    
    return [ Gillespie_one_step, ptmax, pvxp, model_tex, param_tex ]

def run_gillespie( ptrs, nt, nm, pmax, xmax, nthread, nblock, model, param, times, init ):
    #print "\n#### run_gillespie"
    Gillespie_one_step, ptmax, pvxp, model_tex, param_tex = ptrs[0], ptrs[1], ptrs[2], ptrs[3], ptrs[4]
    #print "pointer values, ptmax / pvxp : ", ptmax, pvxp

    # output array
    ntimes = len(times)
    xt_h = np.zeros( [nt,ntimes,xmax] )

    # times array
    t = np.zeros( [nt] ).astype(np.float32)
    d_t = cuda.mem_alloc(t.size*t.dtype.itemsize)
    cuda.memcpy_htod(d_t, t)

    # model 1D array
    d_m = cuda.mem_alloc(model.size*model.dtype.itemsize)
    cuda.memcpy_htod(d_m, model)
    model_tex.set_address(d_m, model.size*model.dtype.itemsize, allow_offset=False)

    # parameters 2D array 
    descr = cuda.ArrayDescriptor()
    descr.width = pmax
    descr.height = nt
    descr.format = cuda.dtype_to_array_format(param.dtype)
    descr.num_channels = 1
    descr.flags = 0

    ary = cuda.Array(descr)

    copy = cuda.Memcpy2D()
    copy.set_src_host(param)
    copy.set_dst_array(ary)
    copy.height = nt
    copy.width_in_bytes = copy.src_pitch = pmax*4
    copy(aligned=True)

    param_tex.set_array(ary)

    # values 2D array, fill with initial conditions
    h_x = np.zeros( [nt,xmax], order="C", dtype=np.float32)
    for i in range(nt):
        for j in range(nm):
            if model[i] == j:
                for k in range(len(init[j])):
                    h_x[i,k] = init[j][k]

    vxw = xmax*4
    alloc_tup = cuda.mem_alloc_pitch(width=vxw, height=nt, access_size=4 )
    vxp = alloc_tup[1]
    d_x = alloc_tup[0]
    #print "\n#### species array pitch :", vxp, "bytes"
    cuda.memcpy_htod(pvxp, np.array([vxp], dtype=np.int32) )

    copy = cuda.Memcpy2D()
    copy.set_src_host(h_x)
    copy.set_dst_device(d_x)
    copy.src_pitch = vxw
    copy.dst_pitch = vxp
    copy.width_in_bytes = vxw
    copy.height = nt
    copy(aligned=True)

    for i in range(nt):
        for k in range(xmax):
            xt_h[i,0,k] = h_x[i,k]

    #print "Running Gillespie"
    for i in range(1,len(times)):
        
        cuda.memcpy_htod(ptmax, ( np.ones( [1] ) * times[i] ).astype(np.float32) )
        #print "\tGillespie step"
        Gillespie_one_step( d_t, d_x, block=(nthread,1,1), grid=(nblock,1) );

        copy = cuda.Memcpy2D()
        copy.set_src_device(d_x)
        copy.set_dst_host(h_x)
        copy.src_pitch = vxp
        copy.dst_pitch = vxw
        copy.width_in_bytes = vxw
        copy.height = nt
        copy(aligned=True)

        #for k in range(nt):
        #    print k, model[k], param[k], h_x[k,0], h_x[k,1]

        for j in range(nt):
            for k in range(xmax):
                xt_h[j,i,k] = h_x[j,k]

    return xt_h

def compile_hazards_test( application_code, mt_data, nthread, nblock ):

    print "\n#### GPU info"
    print "#### max threads / registers / blocks :", tools.DeviceData().max_threads, tools.DeviceData().registers, tools.DeviceData().thread_blocks_per_mp

    f = open(application_code,'r')

    _code_ = _mersenne_twister_source_ + _textures_ + _sample_ + f.read() + _test_hazards_source_
    
    compiled = pycuda.compiler.SourceModule( _code_ )
    pvxp = int(compiled.get_global("vxp")[0])
    
    print "\n#### compile_hazards_test"
    print "#### obtaining vxp pointer of size", compiled.get_global("vxp")[1], "bytes, value :", compiled.get_global("vxp")[0]  

    initialise_twisters( mt_data, compiled, nthread, nblock )

    hazards_test = compiled.get_function("test_hazards")
    model_tex = compiled.get_texref("model_tex")
    param_tex = compiled.get_texref("param_tex")
 
    return [ hazards_test, pvxp, model_tex, param_tex ]

def run_hazards_test( ptrs, nt, nm, pmax, xmax, nthread, nblock, model, param, xin, hout, xout ):
    
    hazards_test, pvxp, model_tex, param_tex = ptrs[0], ptrs[1], ptrs[2], ptrs[3]

    # model 1D array
    d_m = cuda.mem_alloc(model.size*model.dtype.itemsize)
    cuda.memcpy_htod(d_m, model)
    model_tex.set_address(d_m, model.size*model.dtype.itemsize, allow_offset=False)

    # parameters 2D array 
    descr = cuda.ArrayDescriptor()
    descr.width = pmax
    descr.height = nt
    descr.format = cuda.dtype_to_array_format(param.dtype)
    descr.num_channels = 1
    descr.flags = 0

    ary = cuda.Array(descr)

    copy = cuda.Memcpy2D()
    copy.set_src_host(param)
    copy.set_dst_array(ary)
    copy.height = nt
    copy.width_in_bytes = copy.src_pitch = pmax*4
    copy(aligned=True)

    param_tex.set_array(ary)

    # values 2D array, fill with initial conditions
    h_x = np.zeros( [nt,xmax], order="C", dtype=np.float32)
    for i in range(nt):
        for k in range( len(xin) ):
                    h_x[i,k] = xin[k]

    vxw = xmax*4
    alloc_tup = cuda.mem_alloc_pitch(width=vxw, height=nt, access_size=4 )
    vxp = alloc_tup[1]
    d_x = alloc_tup[0]
    cuda.memcpy_htod(pvxp, np.array([vxp], dtype=np.int32) )

    copy = cuda.Memcpy2D()
    copy.set_src_host(h_x)
    copy.set_dst_device(d_x)
    copy.src_pitch = vxw
    copy.dst_pitch = vxp
    copy.width_in_bytes = vxw
    copy.height = nt
    copy(aligned=True)

    # hazards 1D array
    h_h = np.zeros( [nt], dtype=np.int32)
    d_h = cuda.mem_alloc(h_h.size*h_h.dtype.itemsize)
    cuda.memcpy_htod(d_h, h_h)

    hazards_test( d_h, d_x, block=(nthread,1,1), grid=(nblock,1) );

    # copy the x values over
    copy = cuda.Memcpy2D()
    copy.set_src_device(d_x)
    copy.set_dst_host(h_x)
    copy.src_pitch = vxp
    copy.dst_pitch = vxw
    copy.width_in_bytes = vxw
    copy.height = nt
    copy(aligned=True)

    # copy the reactions over
    cuda.memcpy_dtoh(h_h, d_h)
    
    for i in range(nt):
        hout[i] = h_h[i]
        for j in range(xmax):
            xout[i,j] = h_x[i,j]

def compile_cuda_gillespie2( application_code, mt_data, nthread, nblock ):

    print "\n#### GPU info"
    print "#### max threads / registers / blocks :", tools.DeviceData().max_threads, tools.DeviceData().registers, tools.DeviceData().thread_blocks_per_mp

    f = open(application_code,'r')

    _code_ = _mersenne_twister_source_ + _textures_ + _sample_ + f.read() + _gillespie2_source_
    
    compiled = pycuda.compiler.SourceModule( _code_ )
    pvxp = int(compiled.get_global("vxp")[0])

    print "\n#### compile_cuda_gillespie"
    print "#### obtaining vxp pointer of size", compiled.get_global("vxp")[1], "bytes, value :", compiled.get_global("vxp")[0]  

    initialise_twisters( mt_data, compiled, nthread, nblock )

    Gillespie_one_step = compiled.get_function("Gillespie_one_step")
    model_tex = compiled.get_texref("model_tex")
    param_tex = compiled.get_texref("param_tex")

    print "\n#### kernel info" 
    print "#### kernel mem local / shared / registers : ", Gillespie_one_step.local_size_bytes, Gillespie_one_step.shared_size_bytes, Gillespie_one_step.num_regs

    print "\n#### occupancy info"
    occ = tools.OccupancyRecord( tools.DeviceData(), threads=nthread, shared_mem=Gillespie_one_step.shared_size_bytes, registers=Gillespie_one_step.num_regs )
    print "#### threadblocks per mp / limit / occupancy :", occ.tb_per_mp, occ.limited_by, occ.occupancy 
    
    return [ Gillespie_one_step, pvxp, model_tex, param_tex ]

def copy2D_host_to_device(host, dev, src_pitch, dst_pitch, width, height ):
    copy = cuda.Memcpy2D()
    copy.set_src_host(host)
    copy.set_dst_device(dev)
    copy.src_pitch = src_pitch
    copy.dst_pitch = dst_pitch
    copy.width_in_bytes = width
    copy.height = height
    copy(aligned=True)

def copy2D_device_to_host(dev, host, src_pitch, dst_pitch, width, height ):
    copy = cuda.Memcpy2D()
    copy.set_src_device(dev)
    copy.set_dst_host(host)
    copy.src_pitch = src_pitch
    copy.dst_pitch = dst_pitch
    copy.width_in_bytes = width
    copy.height = height
    copy(aligned=True)

def copy2D_host_to_array(host, arr, width, height ):
    copy = cuda.Memcpy2D()
    copy.set_src_host(host)
    copy.set_dst_array(arr)
    copy.height = height
    copy.width_in_bytes = copy.src_pitch = width
    copy.height = height
    copy(aligned=True)

def realize_particles(ptrs, nthread, nblock, xmax, pmax,
                      numOutput, nmodel, model_np, model_nx, dead_models,
                      sampleTheModel, sampleTheParameter, t, priors, kernel, restart, modelDistrib0, nrepeat,
                      euclidianDistance, howToFitData, integrationType, InitValues, timepoints, epsilon, fit, data,
                      margins_prev, model_prev, weights_prev, parameters_prev, 
                      model_curr, parameters_curr, beta_values, hits, sampling):

    nt = nthread*nblock

    Gillespie_one_step, pvxp, model_tex, param_tex = ptrs[0], ptrs[1], ptrs[2], ptrs[3]
    #print "pointer values, ptmax / pvxp : ", ptmax, pvxp

    # output array
    ntimes = len(timepoints)
    xt = np.zeros( [nt,ntimes,xmax] )

    # times array
    h_t = np.zeros( [nt] ).astype(np.float32)
    #h_t = cuda.pagelocked_zeros( [nt], dtype=np.float32)
    d_t = cuda.mem_alloc(h_t.size*h_t.dtype.itemsize)

    # max array
    h_tmax = np.zeros( [nt] ).astype(np.float32)
    #h_tmax = cuda.pagelocked_zeros( [nt], dtype=np.float32)
    d_tmax = cuda.mem_alloc(h_tmax.size*h_tmax.dtype.itemsize)

    # model 1D array
    h_m = np.zeros( [nt], dtype=np.int32 )
    d_m = cuda.mem_alloc(h_m.size*h_m.dtype.itemsize)
   
    # parameters 2D array 
    h_p = np.zeros( [nt,pmax], dtype=np.float32)
    descr = cuda.ArrayDescriptor()
    descr.width = pmax
    descr.height = nt
    descr.format = cuda.dtype_to_array_format(np.float32)
    descr.num_channels = 1
    descr.flags = 0
    d_p = cuda.Array(descr)

    # values 2D matrix
    h_x = np.zeros( [nt,xmax], order="C", dtype=np.float32)
    #h_x = cuda.pagelocked_zeros( [nt,xmax], order="C", dtype=np.float32)
    vxw = xmax*4 #32 bit = 4bytes
    alloc_tup = cuda.mem_alloc_pitch(width=vxw, height=nt, access_size=4 )
    vxp = alloc_tup[1]
    d_x = alloc_tup[0]
    cuda.memcpy_htod(pvxp, np.array([vxp], dtype=np.int32) )

    # Gillespie + distance calculation algorithm
    passed = 0
    failed = 0 
    accepted = 0
    iterations = 0

    attempts = np.zeros([nt], dtype=np.int32)
    status = np.zeros([nt], dtype=np.int32)
    step = np.zeros([nt], dtype=np.int32)
    betas = np.zeros([nt], dtype=np.int32)
    distance = np.zeros([nt], dtype=np.float32)

    while( accepted < numOutput ):
        # refill/reset particles
        # status = 0 : resample particle 
        # status = 1 : reset particle time series
        # status = 2 : particle still valid
        # status = 3 : accepted particle

        #print "iteration:", iterations

        for i in range(nt):
            if(status[i] == 0):
                attempts[i] = 0
                betas[i] = 0
                step[i] = 1
                status[i] = 0
                h_p[i,:] = 0
                h_x[i,:] = 0
                xt[i,:,:] = 0

                m = sampleTheModel(t, nmodel, margins_prev,  modelDistrib0, dead_models, restart )
                p = sampleTheParameter(t, numOutput, priors[m], kernel[m], m, margins_prev, model_prev, weights_prev, parameters_prev, restart ) 

                h_m[i] = m
                for k in range(model_np[m]):
                    h_p[i,k] = p[k]

                #print "sampled h_p:", h_p[i]
                #print "thread:", i, "newparticle"

                for k in range(model_nx[m]):
                    h_x[i,k] = InitValues[m][k]
                    xt[i,0,k] = h_x[i,k]

                h_t[i] = timepoints[0]
                
            if(status[i] == 1):
                step[i] = 1
                status[i] = 0
                h_x[i,:] = 0
                xt[i,:,:] = 0

                for k in range(model_nx[m]):
                    h_x[i,k] = InitValues[m][k]
                    xt[i,0,k] = h_x[i,k]
                
                h_t[i] = timepoints[0]
                
            h_tmax[i] = timepoints[step[i]]
            #print "particle:", h_m[i], h_p[i], h_t[i], h_tmax[i]
 
        cuda.memcpy_htod(d_m, h_m)
        model_tex.set_address(d_m, h_m.size*h_m.dtype.itemsize, allow_offset=False)

        copy2D_host_to_array(h_p, d_p, pmax*4, nt )
        param_tex.set_array(d_p)

        copy2D_host_to_device(h_x, d_x, vxw, vxp, vxw, nt )
        cuda.memcpy_htod(d_t, h_t)
        cuda.memcpy_htod(d_tmax, h_tmax)
        
        #print "Running cuda"
        #print h_t
        #print h_x
        #print h_tmax
        Gillespie_one_step( d_t, d_x, d_tmax, block=(nthread,1,1), grid=(nblock,1) );

        copy2D_device_to_host(d_x, h_x, vxp, vxw, vxw, nt )
        cuda.memcpy_dtoh(h_t, d_t)

        # Do the distance calculation
        for i in range(nt):

            m = h_m[i]
            for k in range(xmax):
                xt[i,step[i],k] = h_x[i,k]
                
            sim = xt[i,0:step[i],0:model_nx[m] ]            
            dat = data[0:step[i]]

            #print "sim:\t", sim.transpose()
            #print "dat:\t", dat.transpose()

            points=howToFitData( fit[m], sim )
            distance[i] = euclidianDistance(points, dat)            
            #print "dist:\t", distance[i]

            # process particle and set the status
            good = False
            if(distance[i] < epsilon[t] and distance[i] >= 0):
                good = True

            end = False
            if( step[i] == ntimes-1 ):
                end = True

            last_attempt = False
            if( attempts[i] == nrepeat-1 ):
                last_attempt = True

            #print "thread:", i, "good / end / last:", good, end, last_attempt

            if(good and not end):
                step[i] = step[i] + 1
                status[i] = 2
                
            elif(good and end and last_attempt):
                attempts[i] = attempts[i] + 1
                betas[i] = betas[i] + 1
                status[i] = 3

            elif(good and end and not last_attempt):
                attempts[i] = attempts[i] + 1
                betas[i] = betas[i] + 1
                status[i] = 1

            elif(not good and last_attempt and betas[i] > 0):
                attempts[i] = attempts[i] + 1
                status[i] = 3

            elif(not good and last_attempt and betas[i] == 0):
                attempts[i] = attempts[i] + 1
                status[i] = 0
                failed = failed + 1

            elif(not good and not last_attempt):
                attempts[i] = attempts[i] + 1
                status[i] = 1

            else:
                print "ABORT: CASE NOT CAUGHT"
                sys.exit()

            #print i, "results pos:", h_tmax[i], h_t[i], h_x[i,:], distance[i], step[i], status[i]
            #print "\n\n"

            # Copy good particles over
            if(status[i] == 3):
                passed = passed + 1

                if(accepted < numOutput):
                    model_curr[accepted] = h_m[i]
                    for k in range(model_np[h_m[i]]):
                        parameters_curr[accepted].append(h_p[i,k])

                    beta_values[accepted] = betas[i]
                    #print "Accepted particle, thread:",i, accepted, model_curr[accepted], parameters_curr[accepted], beta_values[accepted], attempts[i]

                accepted = accepted + 1
                status[i] = 0

        iterations = iterations + 1
    # main loop

    hits[t] = passed
    sampling[t] = passed + failed
    #print "rate:", accepted, sampled, accepted/float(sampled)

    



