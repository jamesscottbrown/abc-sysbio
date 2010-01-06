# cuda utilities

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

    __device__ float rexp(float& lambda, MersenneTwisterState &state, const unsigned int& threadID) {
    return -log( runif(0, 1, state, threadID) )/lambda;
    }
    """

_textures_="""
texture<int, 1, cudaReadModeElementType> model_tex;
texture<float, 2, cudaReadModeElementType> param_tex;
"""

_sample_="""
__device__ int sample(int& nh, float& h0, float* h, MersenneTwisterState &state, const unsigned int& threadID){
int i;
float u = runif(0, h0, state, threadID);
float F = 0;
for(i=0; i<nh; ++i){
F += h[i];
if(F > u) break;
}
return i;
}
"""

_gillespie_source_="""
    //__constant__ float tmax[1];
    //__constant__ int vxp[1];
    __constant__ float tmax;
    __constant__ int vxp;

    __global__ void Gillespie_one_step(float* vt, float* vx) {
    const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

    int i;
    int model = tex1Dfetch(model_tex, tid);

    float x[ NXMAX ];
    //for(i=0; i<nx[model]; ++i){
    for(i=0; i<2; ++i){
       //x[i] = ((float*)( (char*) vx + tid * vxp[0]))[i];
       x[i] = ((float*)( (char*) vx + tid * vxp))[i];
    }
    
    float t  = vt[ tid ];
    float h0;
    
    //while( t < tmax[0] ){
    while( t < tmax ){
       hazards(h0, x, model, MTS[tid], tid);
       if(h0 <= 0.0) break;
       t = t + rexp( h0, MTS[tid], tid );
    }

    vt[tid] = t;
    //for(i=0; i<nx[model]; ++i){
    for(i=0; i<2; ++i){
       //((float*)( (char*) vx + tid * vxp[0]))[i] = x[i];
       ((float*)( (char*) vx + tid * vxp))[i] = x[i];
    }
    }
    """

def initialise_twisters( mt_data, mod, nthread, nblock ):

    print "\n#### initialise_twisters"
    print "#### obtaining MT  pointer of size", mod.get_global("MT")[1], "bytes" 
    print "#### obtaining MTS pointer of size", mod.get_global("MTS")[1], "bytes" 

    pMT = int(mod.get_global("MT")[0])
    
    MT_RNG_COUNT = 32768
    with open(mt_data, 'rb') as f:
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
    

def compile_cuda_gillespie( application_code, mt_data, nthread, nblock ):

    print "\n#### GPU info"
    print "#### max threads / registers / blocks :", tools.DeviceData().max_threads, tools.DeviceData().registers, tools.DeviceData().thread_blocks_per_mp

    f = open(application_code,'r')

    _code_ = _mersenne_twister_source_ + _textures_ + _sample_ + f.read() + _gillespie_source_
    
    compiled = pycuda.compiler.SourceModule( _code_ )
    ptmax = int(compiled.get_global("tmax")[0])
    pvxp = int(compiled.get_global("vxp")[0])

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
    #print "cuda params", type(param)

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
                h_x[i,:] = init[j][:]

    #print np.shape(h_x), h_x
    #print init[0][0], init[0][1]

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

    xt_h[:,0,:] = h_x[:,:]
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

        #print "t:", times[i], "x0 min:", h_x[:,0].min(), "x0 max:", h_x[:,0].max(), "x1 min:", h_x[:,1].min(), "x1 max:", h_x[:,1].max()  

        #for k in range(nt):
        #    print k, model[k], param[k], h_x[k,0], h_x[k,1]

        xt_h[:,i,:] = h_x[:,:]

    return xt_h
