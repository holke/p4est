#@job_name = p4_box_28racks_release_strong_lvl9
#@comment  = "box kegel refined"
#@error    = /homea/hpclab/train056/$(job_name).$(jobid).err
#@output   = /homea/hpclab/train056/$(job_name).$(jobid).out
#@environment = COPY_ALL
#@job_type = bluegene
#@notification = always
#@notify_user = holke@ins.uni-bonn.de
#@wall_clock_limit = 00:15:00
#@bg_size  = 28672
#@queue

bg_size=28672
NUM_RPN=$((16))    # number of processes per node (RPN*THREADS<=64)
NUM_PROC=$(($NUM_RPN*$bg_size)) # number of total processes (bg_size*RPN)
DIR="/work/hbn26/hbn264/exec/p4est/example/ExtremeScaling"
EXEC="$DIR/p8est_bunny"
ARGS="$DIR/p8est_box_tetgen"
echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS

echo "------------------"

bg_size=24576 # 24 racks
NUM_RPN=$((16))    # number of processes per node (RPN*THREADS<=64)
NUM_PROC=$(($NUM_RPN*$bg_size)) # number of total processes (bg_size*RPN)
echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS

echo "------------------"

bg_size=20480 # 20 racks
NUM_RPN=$((16))    # number of processes per node (RPN*THREADS<=64)
NUM_PROC=$(($NUM_RPN*$bg_size)) # number of total processes (bg_size*RPN)
echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS

echo "------------------"
bg_size=16384 # 16 racks
NUM_RPN=$((16))    # number of processes per node (RPN*THREADS<=64)
NUM_PROC=$(($NUM_RPN*$bg_size)) # number of total processes (bg_size*RPN)
echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS

echo "------------------"
