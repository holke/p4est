#@job_name = p4_box_28racks_strong_scale
#@comment  = "box kegel refined"
#@error    = /homea/hpclab/train056/$(job_name).$(jobid).err
#@output   = /homea/hpclab/train056/$(job_name).$(jobid).out
#@environment = COPY_ALL
#@job_type = bluegene
#@notification = always
#@notify_user = holke@ins.uni-bonn.de
#@wall_clock_limit = 00:40:00
#@bg_size  = 28672
#@queue

bg_size=28672
NUM_RPN=$((32))    # number of processes per node (RPN*THREADS<=64)
level=9
# NUM_PROC=$(($NUM_RPN*$bg_size)) # number of total processes (bg_size*RPN)
DIR="/work/hbn26/hbn264/exec/p4est/example/ExtremeScaling"
EXEC="$DIR/p8est_bunny"
ARGS="-f $DIR/p8est_box_tetgen -l$level"

# for i in `seq 16 19`
# do
#   NUM_PROC=$((2**$i)) # number of total processes (bg_size*RPN)
#   echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
#   runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS
#   echo "-------------------------------------------"    
# done

# 8 racks
NUM_PROC=262144
echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS
echo "-------------------------------------------"    

# 16 racks
NUM_PROC=524288
echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS
echo "-------------------------------------------"    


# 20 racks
NUM_PROC=655360 # number of total processes (bg_size*RPN)
echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS
echo "-------------------------------------------"    

# 24 racks
NUM_PROC=786432 # number of total processes (bg_size*RPN)
echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS
echo "-------------------------------------------"    

# 28 racks
NUM_PROC=917504 # number of total processes (bg_size*RPN)
echo "runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS"
runjob --np $NUM_PROC --ranks-per-node=$NUM_RPN : $EXEC $ARGS
echo "-------------------------------------------"    

