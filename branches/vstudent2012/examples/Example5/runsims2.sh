#!/bin/bash

# Need some setup for running cbarnes python
export PATH=/cluster/home/cbarnes/soft/bin/:$PATH
export LD_LIBRARY_PATH=/cluster/home/cbarnes/soft/lib:/cluster/soft/Linux_2.6_64/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:/usr/local/cuda/lib:/cluster/home/cbarnes/soft/lib:${LD_LIBRARY_PATH}

# set python path to read local abcsysbio and cudsim modules
export PYTHONPATH=/cluster/home/filippi/work/Research/Students/MRCProjectABC_JanApril2012/abcsysbio2.06_version03042010_students/abc-sysbio-2.05/
export PYTHONPATH=$PYTHONPATH:/cluster/home/juliane/work/2011/Hess1/cudaSim/cuda-sim-0.06

export PATH=/usr/local/cuda/bin/:/cluster/home/cbarnes/soft/bin:${PATH}


python_exe=/cluster/home/cbarnes/soft/bin/python
abcSysBio_exe=/cluster/home/filippi/work/Research/Students/MRCProjectABC_JanApril2012/abcsysbio2.06_version03042010_students/abc-sysbio-2.05/scripts/run-abc-sysbio


cd /cluster/home/filippi/work/Research/Students/MRCProjectABC_JanApril2012/abcsysbio2.06_version03042010_students/abc-sysbio-2.05/examples/Example5


$python_exe -u $abcSysBio_exe -i Repressilator_submodels_input_file_auto.xml -of=results_5 -cu -lc >& log_all.5

