###################################################
# INSTALLATION OF C LIBRARY OF STATISTICS FUNCTIONS
###################################################
 
1) Update the gsl library paths in the makefile: 

	LIBS = -L/*insert_gsl_library_path

Note that the two libraries to include are gsl and gslcblas 

2) Build the library using the makefile:

	make

3) Move the shared library libstats.so.1.0 into the abcsysbio folder

	mv libstats.co.1.0 ..