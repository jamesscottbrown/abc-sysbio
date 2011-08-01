#!/bin/bash

g++ -c -fPIC newmat1.cpp -o newmat1.o
g++ -c -fPIC newmat2.cpp -o newmat2.o
g++ -c -fPIC newmat3.cpp -o newmat3.o
g++ -c -fPIC newmat4.cpp -o newmat4.o
g++ -c -fPIC newmat5.cpp -o newmat5.o
g++ -c -fPIC newmat6.cpp -o newmat6.o
g++ -c -fPIC newmat7.cpp -o newmat7.o
g++ -c -fPIC newmat8.cpp -o newmat8.o
g++ -c -fPIC newmatex.cpp -o newmatex.o
g++ -c -fPIC bandmat.cpp -o bandmat.o
g++ -c -fPIC submat.cpp -o submat.o
g++ -c -fPIC myexcept.cpp -o myexcept.o
g++ -c -fPIC cholesky.cpp -o cholesky.o
g++ -c -fPIC evalue.cpp -o evalue.o
g++ -c -fPIC fft.cpp -o fft.o
g++ -c -fPIC hholder.cpp -o hholder.o
g++ -c -fPIC jacobi.cpp -o jacobi.o
g++ -c -fPIC newfft.cpp -o newfft.o
g++ -c -fPIC sort.cpp -o sort.o
g++ -c -fPIC svd.cpp -o svd.o
g++ -c -fPIC nm_misc.cpp -o nm_misc.o
g++ -c -fPIC newmatrm.cpp -o newmatrm.o
g++ -c -fPIC newmat9.cpp -o newmat9.o

platform=`uname`

if [[ "$platform" == "Darwin" ]]
then
    loc=`pwd`
    g++ -dynamiclib *.o -o libnewmat.so;
    install_name_tool -id ${loc}/libnewmat.so ${loc}/libnewmat.so
else
    g++ -shared *.o -o libnewmat.so
fi
