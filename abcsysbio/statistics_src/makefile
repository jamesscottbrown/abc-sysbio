
INCS = -I/cluster/soft/Linux_2.6_64_project/include
LIBS = -L/cluster/soft/Linux_2.6_64_project/lib -lgsl -lgslcblas
CXXFLAGS = -fPIC

prog:
	g++ $(CXXFLAGS) -c statistics.cpp -o statistics.o $(INCS)
	g++ -shared -Wl -o libstats.so.1.0 statistics.o $(LIBS)

clean:
	rm -f *.o



