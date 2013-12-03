# SET UP YOUR OWN PATHS HERE
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/dev/abc-sysbio-area/abc-sysbio/branches/structure/
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/dev/cuda-sim-area/cuda-sim/trunk/

exe=/home/cbarnes/dev/abc-sysbio-area/abc-sysbio/branches/structure/scripts/run-abc-sysbio

# simulate
#python -u ${exe} -i input_asym2_rep_gen1_sde.xml -f -of=run_sim --cuda --localcode -S --timing --custd=dist_pulseB-v1

# fit 
python -u ${exe} -i input_asym2_rep_gen1_sde.xml -f -of=run_fit --cuda --localcode --design --timing --custd=dist_pulseB-v1
