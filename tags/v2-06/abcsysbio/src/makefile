
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
LIBARG := -shared
TARG := -lrt
GCC :=
endif
ifeq ($(UNAME), Darwin)
LIBARG := -dynamiclib
TARG :=
GCC := -L/usr/lib -lgcc
endif

LIBS1= -L$(ABC_GSL_LIB) -lgsl -lgslcblas $(TARG)
LIBS2= -L$(ABC_NM_LIB) -lnewmat 
LIBS= $(GCC) $(LIBS1) $(LIBS2)
LOCS= -Wl,-rpath,$(ABC_GSL_LIB) -Wl,-rpath,$(ABC_NM_LIB) -Wl,-rpath,. 

INCLUDES= -I. -I$(ABC_SRC_DIR) -I$(ABC_GSL_INC) -I$(ABC_NM_INC)

CXXFLAGS= -fPIC -w -O2

#all :  model solver_base sde_solver_base library cleanup
all : prep library cleanup

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
	g++ $(CXXFLAGS) $(LIBARG) -o $(LIBNAME) ChildModel.o Solver.o SDESolver.o $(SOLVER).o $(LIBS) $(LOCS) 

cleanup :
	rm ChildModel.cpp ChildModel.hpp ChildModel.o Solver.o SDESolver.o $(SOLVER).o

prep :
	rm -f $(LIBNAME)