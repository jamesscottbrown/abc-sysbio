# testing python package cudatools

import numpy
from numpy.random import *
import abccuda
import abccuda.gillespie, abccuda.lsoda
import gillespie, lsoda
import time, csv, string

# get the location of the MersenneTwister.cu file
s = abccuda.__file__
mt_cu = s.replace('__init__.pyc','MersenneTwister.cu')
print '#### using MersenneTwister.cu at', mt_cu

# get the location of the MersenneTwister.dat file
s = abccuda.__file__
mt_data = s.replace('__init__.pyc','MersenneTwister.dat')
print '#### using MersenneTwister.dat at', mt_data

# get the location of the cuLsoda_all.cu file
s = abccuda.__file__
ls_cu = s.replace('__init__.pyc','cuLsoda_all.cu')
print '#### using cuLsoda_all.cu at', ls_cu

# define some printing functions
def print_mat(u,nt,nb,nr):
    for i in range(nt*nb):
        for j in range(nr):
            print u[j + i*nr],
        print ''

def print_gill_results(sims,nt,xmax,times,outfile):
    f = open(outfile,'w')
    for k in range(nt):
        for kk in range(xmax):
            print >>f, kk,
            for kkk in range(len(times)):
                print >>f, ",", sims[k,kkk,kk].astype(numpy.int32), 
            print >>f, "\n",  
    f.close()

def print_ode_results(sims,nt,xmax,times,outfile):
    f = open(outfile,'w')
    for k in range(nt):
        for kk in range(xmax):
            print >>f, kk,
            for kkk in range(len(times)):
                print >>f, ",", sims[k,kkk,kk].astype(numpy.float64), 
            print >>f, "\n",  
    f.close()

if(1):
    nthread = 32
    nblock = 8
    nr = 10

    # test the persistence of the MersenneTwister
    numpy.random.seed(10) 
    c = gillespie.compile_mt_code( gillespie._rng_test_source_, mt_cu, mt_data, nthread, nblock, options=['--maxrregcount=30'] )
    print '\n\n\n##### run1 #####'
    u = gillespie.run_mt_test(c, nthread, nblock, nr )  
    print_mat( u, nthread, nblock, nr )
    print '\n\n\n##### run2 #####'
    u = gillespie.run_mt_test(c, nthread, nblock, nr )  
    print_mat( u, nthread, nblock, nr )

    # recompile and run to make sure we get the same numbers
    numpy.random.seed(10)
    c = gillespie.compile_mt_code( gillespie._rng_test_source_, mt_cu, mt_data, nthread, nblock, options=['--maxrregcount=30'] )
    print '\n\n\n##### run12 #####'
    u = gillespie.run_mt_test(c, nthread, nblock, nr*2 )  
    print_mat( u, nthread, nblock, nr*2 )


if(1):
    numpy.random.seed(9)

    nthread = 32
    nblock = 8

    appfile = "../examples/Example4/immdeath_mjp.cu"
    times = numpy.arange(0,21)
    init = [[1]]

    nm = 1
    xmax = 1
    pmax = 2
    
    c = gillespie.compile_gillespie(appfile, mt_cu, mt_data, nthread, nblock, pmax, xmax, options=[], write=False )

    model = numpy.array( [0 for i in range(nthread*nblock)], dtype=numpy.int32)
    param = numpy.zeros( [nthread*nblock, 2], dtype=numpy.float32)
    param[:,0] = 1.0
    param[:,1] = 0.1
   
    start_time = time.time()
    u = gillespie.run_gillespie(c, nthread, nblock, nm, pmax, xmax,  model, param, times, init )
    print '#### timing:', time.time() - start_time

    print_gill_results(u,nthread*nblock,xmax,times,'immdeath_mjp.tout')

if(1):
    numpy.random.seed(9)
    
    nthread = 32
    nblock = 8

    appfile = "../examples/Example4/immdeath_ode.cu"
    times = numpy.arange(0,21)
    init = [[1]]

    nm = 1
    xmax = 1
    pmax = 2
    dt = -1
    
    #options=["-Xopencc","-O0"]
    c = lsoda.compile_lsoda( appfile, ls_cu, nthread, nblock, pmax, xmax, options=[], write=False )

    model = numpy.array( [0 for i in range(nthread*nblock)], dtype=numpy.int32)
    model_nx = [1]

    param = numpy.zeros( [nthread*nblock, pmax], dtype=numpy.float32)
    param[:,0] = 1.0
    param[:,1] = 0.1
     
    start_time = time.time()
    u = lsoda.run_lsoda(c, nthread, nblock, nm, pmax, xmax,  model_nx, model, param, times, init, dt )
    print '#### timing:', time.time() - start_time

    print_ode_results(u,nthread*nblock,xmax,times,'immdeath_ode.tout')

if(1):
    numpy.random.seed(9)
    
    nthread = 32
    nblock = 8
    pmax = 6
    xmax = 4
    nm = 3
    dt = -1
    init = [ [20.0, 10.0, 0.0 ], [20.0, 0.0, 10.0, 0.0 ], [20.0, 10.0, 0.0] ]
    
    appfile = '../examples/Example1/SIRmodels.cu'

    #options=["-Xopencc","-O0"]
    c = lsoda.compile_lsoda( appfile, ls_cu, nthread, nblock, pmax, xmax, options=[], write=False )

    times = numpy.array( [ 0.0, 0.6, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 ] )
                      
    model = numpy.zeros( [nthread*nblock], dtype=numpy.int32)
    param = numpy.zeros( [nthread*nblock, pmax], dtype=numpy.float32)

    model_nx = [3,4,3]

    for i in range(nthread*nblock):
        model[i] = numpy.random.randint(0,3)

        param[i,0] = 1.0
        param[i,1] = numpy.random.uniform(0,5)
        param[i,2] = numpy.random.uniform(0,5)
        param[i,3] = numpy.random.uniform(0,5)
        param[i,4] = numpy.random.uniform(0,5)
        param[i,5] = numpy.random.uniform(0,10)
      
    start_time = time.time()
    u = lsoda.run_lsoda(c, nthread, nblock, nm, pmax, xmax, model_nx, model, param, times, init, dt )
    print '#### timing:', time.time() - start_time

    print_ode_results(u,nthread*nblock,xmax,times,'SIRmodels_ode.tout')
