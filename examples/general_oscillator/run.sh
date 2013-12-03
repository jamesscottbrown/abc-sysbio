# SET UP YOUR OWN PATHS HERE
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/dev/abc-sysbio-area/abc-sysbio/branches/structure/
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/dev/cuda-sim-area/cuda-sim/trunk/

exe=/home/cbarnes/dev/abc-sysbio-area/abc-sysbio/branches/structure/scripts/run-abc-sysbio

# three node general network
python -u ${exe} --design -i input_three_node_gen_sde_1.xml -f -of=run_gen_three --cuda --localcode --design --timing --custd=dist_pulseB-v1

# three node general network with looser priors
# python -u ${exe} --design -i input_three_node_gen_sde_1_loose.xml -f -of=run_gen_three_loose --cuda --localcode --design --timing --custd=dist_pulseB-v1

# two node network
# python -u ${exe} --design -i input_two_node_gen_sde_1.xml -f -of=run_gen_two --cuda --localcode --design --timing --custd=dist_pulseB-v1