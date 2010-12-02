# This makefile can be used to run all the examples
# Edit the variables to suit your needs

# VERSION OF THE PACKAGE TO TEST
VERSION := 1.02

# LOCATION OF THE DOWNLOADED TAR BALL
SOURCE := ../dev/abc-sysbio-area/abc-sysbio/dist/abc-sysbio-$(VERSION).tar.gz

# LOCATION OF RUN_ABC_SYSBIO
RUN_ABC := bin/run-abc-sysbio

# THE TEST DIRECTORY TO INSTALL INTO
# ASSUMES YOU ARE IN THIS DIRECTORY
TEST := /cluster/home/cbarnes/test

all: install setup ex2 ex4 ex1 ex3

install:
	cp $(SOURCE) .
	tar -xzf abc-sysbio-$(VERSION).tar.gz
	cd abc-sysbio-$(VERSION); python setup.py install --prefix=$(TEST)
	export PYTHONPATH=lib/python2.6/site-packages/

setup:	
	cp abc-sysbio-$(VERSION)/examples/Example1/* .
	cp abc-sysbio-$(VERSION)/examples/Example2/* .
	cp abc-sysbio-$(VERSION)/examples/Example4/* .	
	cp abc-sysbio-$(VERSION)/examples/Example3/* .

ex1:
	python -u $(RUN_ABC) --infile input_file_SIR.txt -f -of=res_ex1_rm0 --timing -rm=0 > log.ex1_rm0
	python -u $(RUN_ABC) --infile input_file_SIR.txt -f -of=res_ex1_rm1 --timing -rm=1 > log.ex1_rm1
	python -u $(RUN_ABC) --infile input_file_SIR.txt -f -of=res_ex1_rm2 -rm=2 -ct=64 -cb=32 --cudacode=SIRmodels.cu -pm --timing > log.ex1_rm1

ex2:
	python -u $(RUN_ABC) --infile input_file_lingrow.txt -f -of=res_ex2_rm1 --timing -rm=1 -pm > log.ex2_rm1
	python -u $(RUN_ABC) --infile input_file_lingrow.txt -f -of=res_ex2_rm0 --timing -rm=0 -pm > log.ex2_rm0

ex4:
	python -u $(RUN_ABC) -sd=10 --infile input_file_pop.txt -sd=10 -f --runmode=2 -of=res_ex4_rm2 -ct=256 -cb=30 --cudacode=immdeath_logistic.cu -pm -t > log.ex4_rm2

ex3:	
	python -u $(RUN_ABC) --infile input_file_abcRejection.txt -f -of=res_ex3_rm1 --timing -rm=1 > log.ex3_rm1


clean:
	rm -f xml *.pyc *.py input_file_SIR.txt _data.png
