##export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/soft/Linux_2.6_64/lib:/cluster/home/cbarnes/soft/code/newmat11/:.

LIBS1= -L$(ABC_GSL_LIB) -lgsl -lgslcblas -lrt
LIBS2= -L$(ABC_NM_LIB) -lnewmat
LIBS= $(LIBS1) $(LIBS2)
INCLUDES= -I. -I$(ABC_SRC_DIR) -I$(ABC_GSL_INC) -I$(ABC_NM_INC)

CXXFLAGS= -fPIC -w

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

solver_base : $(ABC_SRC_DIR)/Solver.cpp
	g++ $(CXXFLAGS) -c $(ABC_SRC_DIR)/Solver.cpp -o Solver.o $(INCLUDES)

sde_solver_base : solver_base $(ABC_SRC_DIR)/SDESolver.cpp
	g++ $(CXXFLAGS) -c $(ABC_SRC_DIR)/SDESolver.cpp -o SDESolver.o $(INCLUDES)

solver : solver_base
	g++ $(CXXFLAGS) -c $(ABC_SRC_DIR)/$(SOLVER).cpp -o $(SOLVER).o $(INCLUDES)

library: model solver_base sde_solver_base solver
	g++ $(CXXFLAGS) -shared -Wl,-soname,$(LIBNAME) -o $(LIBNAME) ChildModel.o Solver.o SDESolver.o $(SOLVER).o $(LIBS)

cleanup :
	rm ChildModel.cpp ChildModel.hpp ChildModel.o Solver.o $(SOLVER).o $(LIBNAME)

