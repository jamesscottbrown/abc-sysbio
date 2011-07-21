
##MODEL=p53IVModel
##SOLVER=EulerSDESolver
##LIBNAME=$(MODEL).so.1.0

LIBS1= -L/project/soft/linux/lib/ -lgsl -lgslcblas -lrt
LIBS2= -L/project/home/student3/Desktop/newmat11 -lnewmat
INCLUDES1= -I/project/soft/linux/include/ 
INCLUDES2= -I/project/home/student3/Desktop/newmat11 -I.
CXXFLAGS= -fPIC

all : clean debug model_base model solver_base sde_solver_base   solver_plus_model library #include sde_solver_base

debug:
	echo "\n\n\n\n\ MAKEFILE"
	echo $(MODEL)
	echo $(SOLVER)
	echo $(LIBNAME)
	echo $(CHILDMODEL)

model_base : Model.cpp Model.hpp	
	g++ $(CXXFLAGS) -c Model.cpp -o Model.o $(INCLUDES1) $(INCLUDES2)

model : model_base
	cp $(MODEL).cpp ChildModel.cpp
	cp $(MODEL).hpp ChildModel.hpp
	g++ $(CXXFLAGS) -c ChildModel.cpp -o ChildModel.o $(INCLUDES1) $(INCLUDES2)

solver_base : Solver.cpp Solver.hpp
	g++ $(CXXFLAGS) -c Solver.cpp -o Solver.o $(INCLUDES1) $(INCLUDES2)

sde_solver_base : SDESolver.cpp SDESolver.hpp
	g++ $(CXXFLAGS) -c SDESolver.cpp -o SDESolver.o $(INCLUDES1) $(INCLUDES2) 


solver_plus_model : model_base solver_base sde_solver_base model # include sde_solver_base
	g++ $(CXXFLAGS) -c $(SOLVER).cpp -o $(MODEL).o SDESolver.o $(INCLUDES1) $(INCLUDES2)

library: model_base model solver_base sde_solver_base  solver_plus_model #include sde_solver_base
	g++ $(CXXFLAGS) -shared -Wl,-soname,$(LIBNAME) -o $(LIBNAME) Model.o ChildModel.o SDESolver.o Solver.o $(MODEL).o $(LIBS1) $(LIBS2)
##include SDESolver.o if you want it to run for SDE algorithms
##model_base solver_base model solver_plus_model
##Model.o Solver.o ChildModel.o $(MODEL).o

clean :
	rm -f *.o lib*


