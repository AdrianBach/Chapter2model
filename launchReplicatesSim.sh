#!/bin/usr/env bash

# Author : Adrian Bach
# Date : 4th Feb 2022

# To use this correctly, the directory must contain:
# - this script
# - the latest version of the ibm
# - the makefile

# to eventually get access 
#   chmod u+x ./lauchreplicatesSim.sh
# and execute
#   bash ./lauchreplicatesSim.sh


#### Parameter list: assign values ####

## User defined variables ##

# landscape variables
size=5       # argv[2] world's side size
res_nb=2     # argv[3] number of resource types
max_res_1=10 # argv[4] max resource 1 per cell
max_res_2=10 # argv[5] max resource 2 per cell

max_cell=2.5 # max expected number of preys of each kind per cell # if 2.5, between 2 and 3 animals per cell

# prey variables
pry_nb=2        # argv[6] number of prey types
pry_init_1=20   # argv[7] prey 1 initial density in nb of individuals
pry_init_2=20   # argv[8] prey 2 initial density
pry_move_1=0.1  # argv[9] prey 1 max movement range in fraction of size
pry_move_2=0.1  # argv[10] prey 2 movement range
# pry_surv_1=3  # argv[13] prey 1 resource units needed to pass survival trial
# pry_surv_2=3  # argv[14] prey 2 resource units needed to pass survival trial
pry_offs_1=2    # argv[15] prey 1 max number of offspring
pry_offs_2=2    # argv[16] prey 2 max number of offspring
# pry_repr_1=5  # argv[17] prey 1 resource units needed to pass reproduction trial
# pry_repr_2=5  # argv[18] prey 2 resource units needed to pass reproduction trial

# predator variables
prd_nb=1        # argv[19] number of predator types
prd_init_1=5    # argv[20] predator 1 initial density in nb of individuals
prd_move_1=0.1  # argv[21] predator 1 max movement range in fraction of size
prd_offs_1=1    # argv[24] predator 1 max number of offspring
prd_intr_1=0    # argv[26] predator 1 time of introduction in the model
prd_asym_1=0.5  # argv[27] asymmetry in prey1 to prey2 conversion rates
prd_ctch_1=0.5  # argv[28] predator catch probability
prd_oprt_1=1    # argv[29] is predator oportunistic? (0 or 1)

# time variables
simu_time=20    # argv[30] simulation time
freq_repr=10    # argv[31] frequency of reproduction trials
freq_surv=$freq_repr    # argv[32] frequency of survival trials
freq_rfll=$freq_repr    # argv[33] frequency of landscape resources refill

# frequency of assessment
freq_rslt=1    # argv[34] frequency of landscape results shot
freq_snap=101  # argv[35] frequency of snap measure

# number of replicates
rep=1

## Non user defined variables ##

# pry_cons_1= # argv[11] prey 1 max consumption in resource units
divide=$max_res_1; by=$max_cell; pry_cons_1=`echo "scale=0; ($divide+$by-1)/$by" | bc`; # echo "pry_cons_1 = $pry_cons_1" # bash way to ceil a float. by = max expected number of animal per cell HARD CODED NOT IDEAL.
# pry_cons_2= # argv[12] prey 2 max consumption
divide=$max_res_2; by=$max_cell; pry_cons_2=`echo "scale=0; ($divide+$by-1)/$by" | bc`; # echo "pry_cons_2 = $pry_cons_2"
# pry_surv_1= # argv[13] prey 1 resource units needed to pass survival trial. 
divide=$freq_surv*$pry_cons_1; by=3; pry_surv_1=`echo "scale=0; ($divide+$by-1)/$by" | bc`; # echo "pry_surv_1 = $pry_surv_1" # by = a third of the max number of consecutive fasting days.
# pry_surv_2= # argv[14] prey 2 resource units needed to pass survival trial
divide=$freq_surv*$pry_cons_2; by=3; pry_surv_2=`echo "scale=0; ($divide+$by-1)/$by" | bc`; # echo "pry_surv_2 = $pry_surv_2"
pry_repr_1=$pry_surv_1; # argv[17] prey 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.
pry_repr_2=$pry_surv_2; # argv[18] prey 2 resource units needed to pass reproduction trial

prd_cons_1=$((3*$pry_cons_1)) # arg[22]
divide=$prd_cons_1*$freq_surv; by=3; prd_surv_1=`echo "scale=0; ($divide+$by-1)/$by" | bc`; # arg[23]
prd_repr_1=$prd_surv_1; # echo "prd_repr_1 = $prd_repr_1" # argv[25] predator 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.

# name the simulation with only the variables of interest and their value
sim_name="test-size$size-simTime$simu_time-res1max$max_res_1-res2max$max_res_2-pry1init$pry_init_1-pry1cons$pry_cons_1-prdInit$prd_init_1-prdSurv$prd_surv_1-prdCtch$prd_ctch_1" # argv[1]


#### Simulation loop ####

echo "prey populations initial density are $pry_init_1 and $pry_init_2 and predators $prd_init_1"

# Use make to compile the ibm
make

## Create useful directories ##

# if the directory exists, delete content? otherwise create it
if [ -d "$sim_name" ]
then
    while true; do
        read -p "The simulation directory already exists, its content will be deleted, continue?" yn
        case $yn in
            [Yy]* ) find $sim_name -mindepth 1 -delete; mkdir $sim_name; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
else
    # create a directory for the simulation
    mkdir $sim_name
fi

# move the executable program there
cp test-chapter2ibm.o $sim_name/test-chapter2ibm.o

# move to this directory
cd $sim_name

# write a txt file with all the parameter values
printf "sim_name = $sim_name \t # argv[1] \n\n" >> paramFile.txt
printf "# landscape variables\n" >> paramFile.txt
printf "size = $size \t # argv[2] world's side size\n" >> paramFile.txt
printf "res_nb = $res_nb \t # argv[3] number of resource types\n" >> paramFile.txt
printf "max_res_1 = $max_res_1 \t # argv[4] max resource 1 per cell\n" >> paramFile.txt
printf "max_res_2 = $max_res_2 \t # argv[5] max resource 1 per cell\n\n" >> paramFile.txt
printf "# prey variables\n" >> paramFile.txt
printf "pry_nb = $pry_nb \t # argv[6] number of prey types\n" >> paramFile.txt
printf "pry_init_1 = $pry_init_1 \t # argv[7] prey 1 initial density in nb of individuals\n" >> paramFile.txt
printf "pry_init_2 = $pry_init_2 \t # argv[8] prey 2 initial density\n" >> paramFile.txt
printf "pry_move_1 = $pry_move_1 \t # argv[9] prey 1 max movement range in fraction of size\n" >> paramFile.txt
printf "pry_move_2 = $pry_move_2 \t # argv[10] prey 2 movement range\n" >> paramFile.txt
printf "pry_cons_1 = $pry_cons_1 \t # argv[11] prey 1 max consumption in resource units\n" >> paramFile.txt
printf "pry_cons_2 = $pry_cons_2 \t # argv[12] prey 2 max consumption\n" >> paramFile.txt
printf "pry_surv_1 = $pry_surv_1 \t # argv[13] prey 1 resource units needed to pass survival trial\n" >> paramFile.txt
printf "pry_surv_2 = $pry_surv_2 \t # argv[14] prey 2 resource units needed to pass survival trial\n" >> paramFile.txt
printf "pry_offs_1 = $pry_offs_1 \t # argv[15] prey 1 max number of offspring\n" >> paramFile.txt
printf "pry_offs_2 = $pry_offs_2 \t # argv[16] prey 2 max number of offspring\n" >> paramFile.txt
printf "pry_repr_1 = $pry_repr_1 \t # argv[17] prey 1 resource units needed to pass reproduction trial\n" >> paramFile.txt
printf "pry_repr_2 = $pry_repr_2 \t # argv[18] prey 2 resource units needed to pass reproduction trial\n\n" >> paramFile.txt
printf "# predator variables\n" >> paramFile.txt
printf "prd_nb = $prd_nb \t # argv[19] number of predator types\n" >> paramFile.txt
printf "prd_init_1 = $prd_init_1 \t # argv[20] predator 1 initial density in nb of individuals\n" >> paramFile.txt
printf "prd_move_1 = $prd_move_1 \t # argv[21] predator 1 max movement range in fraction of size\n" >> paramFile.txt
printf "prd_cons_1 = $prd_cons_1 \t # argv[22] predator 1 max movement range in fraction of size\n" >> paramFile.txt
printf "prd_surv_1 = $prd_surv_1 \t # argv[23] predator 1 resource units needed to pass survival trial\n" >> paramFile.txt
printf "prd_offs_1 = $prd_offs_1 \t # argv[24] predator 1 max number of offspring\n" >> paramFile.txt
printf "prd_repr_1 = $prd_repr_1 \t # argv[25] predator 1 resource units needed to pass reproduction trial\n" >> paramFile.txt
printf "prd_intr_1 = $prd_intr_1 \t # argv[26] predator 1 time of introduction in the model\n" >> paramFile.txt
printf "prd_asym_1 = $prd_asym_1 \t # argv[27] predator 1 asymmetry in prey1 to prey2 conversion rates\n" >> paramFile.txt
printf "prd_ctch_1 = $prd_ctch_1 \t # argv[28] predator 1 predator catch probablility \n" >> paramFile.txt
printf "prd_gnrl_1 = $prd_oprt_1 \t # argv[29] predator 1 oportunistic? (0 or 1) \n\n" >> paramFile.txt
printf "# time variables\n" >> paramFile.txt
printf "simu_time = $simu_time \t # argv[30] simulation time\n" >> paramFile.txt
printf "freq_repr = $freq_repr \t # argv[31] frequency of reproduction trials\n" >> paramFile.txt
printf "freq_surv = $freq_surv \t # argv[32] frequency of survival trials\n" >> paramFile.txt
printf "freq_rfll = $freq_rfll \t # argv[33] frequency of results measure\n\n" >> paramFile.txt
printf "# frequency of assessment\n" >> paramFile.txt
printf "freq_rslt = $freq_rslt \t # argv[34] frequency of landscape snap shot\n" >> paramFile.txt
printf "freq_snap = $freq_snap \t # argv[35] frequency of results measure\n\n" >> paramFile.txt
printf "# number of replicates\n" >> paramFile.txt
printf "rep = $rep\n\n" >> paramFile.txt
printf "# Simulation infos \n\n" >> paramFile.txt
printf "rep \t seed (arg[36]) \t sim time (s) \t sim time (h) \n\n" >> paramFile.txt

## Start simulation replicates ##

# launch replicates iteratively
for ((j=0 ; j<$rep ; j++))
do
    echo "lauching replicate #$j"

    # seed for random number generator
    rand_seed=$RANDOM # argv[36] set the seed with a random number 
    # rand_seed=$(date +%s) # argv[36] set the seed with a random number 

    start=$(date +%s)
    ./test-chapter2ibm.o $sim_name $size $res_nb $max_res_1 $max_res_2 $pry_nb $pry_init_1 $pry_init_2 $pry_move_1 $pry_move_2 $pry_cons_1 $pry_cons_2 $pry_surv_1 $pry_surv_2 $pry_offs_1 $pry_offs_2 $pry_repr_1 $pry_repr_2 $prd_nb $prd_init_1 $prd_move_1 $prd_cons_1 $prd_surv_1 $prd_offs_1 $prd_repr_1 $prd_intr_1 $prd_asym_1 $prd_ctch_1 $prd_oprt_1 $simu_time $freq_repr $freq_surv $freq_rfll $freq_rslt $freq_snap $rand_seed
    end=$(date +%s)

    time_s=$(($end-$start))
    time_h=$(($time_s/3600))

    printf "$j \t $rand_seed \t\t $time_s \t\t $time_h \n" >> paramFile.txt
    
    mv $sim_name-ResultsTable.csv rep$j-$sim_name-ResultsTable.csv
    mv $sim_name-SnapshotTable.csv rep$j-$sim_name-SnapshotTable.csv
done

# no need for the executable afterwards
rm ./test-chapter2ibm.o

echo "All replicates done, moving back to the upper directory"

# move back to the upper directory
cd ..