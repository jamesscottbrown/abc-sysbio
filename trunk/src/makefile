##export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/soft/Linux_2.6_64/lib:/cluster/home/cbarnes/soft/code/newmat11/:.

LIBS1= -L/cluster/soft/Linux_2.6_64/lib -lgsl -lgslcblas -lrt
LIBS2= -L/cluster/home/cbarnes/soft/code/newmat11/ -lnewmat
LIBS= $(LIBS1) $(LIBS2)

DIR=/cluster/home/cbarnes/dev/abc-sysbio-area/abc-sysbio/trunk/src
INCLUDES= -I/project/soft/linux/include/ -I/cluster/soft/Linux_2.6_64/include -I/cluster/home/cbarnes/soft/code/newmat11/ -I. -I$(DIR)

CXXFLAGS= -fPIC

#all :  model solver_base sde_solver_base library cleanup
all : library

debug:
	echo "MAKEFILE DEBUG"
	echo "MODEL" $(MODEL)
	echo "SOLVER" $(SOLVER)
	echo "LIBNAME" $(LIBNAME)

model : $(MODEL).cpp $(MODEL).hpp
	cp $(MODEL).cpp ChildModel.cpp
	cp $(MODEL).hpp ChildModel.hpp
	g++ $(CXXFLAGS) -c ChildModel.cpp -o ChildModel.o $(INCLUDES)

solver_base : $(DIR)/Solver.cpp
	g++ $(CXXFLAGS) -c $(DIR)/Solver.cpp -o Solver.o $(INCLUDES)

sde_solver_base : solver_base $(DIR)/SDESolver.cpp
	g++ $(CXXFLAGS) -c $(DIR)/SDESolver.cpp -o SDESolver.o $(INCLUDES)

solver : solver_base
	g++ $(CXXFLAGS) -c $(DIR)/$(SOLVER).cpp -o $(SOLVER).o $(INCLUDES)

library: model solver_base sde_solver_base solver
	g++ $(CXXFLAGS) -shared -Wl,-soname,$(LIBNAME) -o $(LIBNAME) ChildModel.o Solver.o SDESolver.o $(SOLVER).o $(LIBS)

cleanup :
	rm ChildModel.cpp ChildModel.hpp ChildModel.o Solver.o $(SOLVER).o

