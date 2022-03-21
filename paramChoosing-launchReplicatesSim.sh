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

#### Launch simulations ####

# Use make to compile the ibm
make

#### Parameter list: assign values ####

## User defined variables ##

# landscape variables
size=4 # argv[2] world's side size
res_nb=2 # argv[3] number of resource types
max_res_1=10 # argv[4] max resource 1 per cell
max_res_2=10 # argv[5] max resource 1 per cell

# prey variables
pry_nb=2 # argv[6] number of prey types
pry_move_1=0.25 # argv[9] prey 1 max movement range in fraction of size
pry_move_2=0.25 # argv[10] prey 2 movement range
pry_offs_1=1 # argv[15] prey 1 max number of offspring
pry_offs_2=1 # argv[16] prey 2 max number of offspring

# predator variables
prd_nb=1 # argv[19] number of predator types
prd_move_1=0.25 # argv[21] predator 1 max movement range in fraction of size
prd_cons_1=1 # argv[22] predator 1 max catches per time step in prey unit
prd_offs_1=1 # argv[25] predator 1 max number of offspring
prd_intr_1=99 # argv[27] predator 1 time of introduction in the model

# time variables
simu_time=10 # argv[28] simulation time
burn_time=0 # argv[29] time steps to burn before lauching regular survival and reproduction trials
freq_repr=2 # argv[30] frequency of reproduction trials
freq_surv=2 # argv[31] frequency of survival trials

# frequency of assessment
freq_rslt=1 # argv[32] frequency of landscape snap shot
freq_snap=5 # argv[33] frequency of results measure

# number of replicates
rep=3

## Non user defined variables ##

# pry_cons_1= # argv[11] prey 1 max consumption in resource units
divide=$max_res_1; by=3; ((pry_cons_1=($divide+$by-1)/$by)); #echo "pry_cons_1 = $pry_cons_1" # bash way to ceil a float. by = max expected number of animal per cell.
# pry_cons_2= # argv[12] prey 2 max consumption
divide=$max_res_2; by=3; ((pry_cons_2=($divide+$by-1)/$by)); # echo "pry_cons_2 = $pry_cons_2"
# pry_surv_1= # argv[13] prey 1 resource units needed to pass survival trial. 
divide=$freq_repr*$pry_cons_1; by=3; ((pry_surv_1=($divide+$by-1)/$by)); # echo "pry_surv_1 = $pry_surv_1" # by = max number of time step without eating.
# pry_surv_2= # argv[14] prey 2 resource units needed to pass survival trial
divide=$freq_repr*$pry_cons_2; by=3; ((pry_surv_2=($divide+$by-1)/$by)); # echo "pry_surv_2 = $pry_surv_2"
pry_repr_1=$pry_surv_1; # echo "pry_repr_1 = $pry_repr_1" # argv[17] prey 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.
pry_repr_2=$pry_surv_2; # echo "pry_repr_2 = $pry_repr_2" # argv[18] prey 2 resource units needed to pass reproduction trial

# prd_conv_1= # argv[23] predator 1 conversion rate of 1 catch into resources
divide=$pry_surv_1*3; by=$freq_surv*$prd_cons_1; ((prd_conv_1=($divide+$by-1)/$by)); # echo "prd_conv_1 = $prd_conv_1" # defined such that prd_surv is not too far from pry_surv_1.
# prd_surv_1= # argv[24] predator 1 resource units needed to pass survival trial
divide=$freq_repr*$prd_cons_1*$prd_conv_1; by=3; ((prd_surv_1=($divide+$by-1)/$by)); # echo "prd_surv_1 = $prd_surv_1" # by = max number of time step without eating.
prd_repr_1=$prd_surv_1; # echo "prd_repr_1 = $prd_repr_1" # argv[26] predator 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.

#### Simulation loop ####

# max=$size*$size*3*2 # 2 times the max expected number of animal per cell times the number of cells
# min=$size*$size*3/10 # one tenth of the maximum

pop_array=($(seq 10 10 40)) # hard coded for now

# echo "pop_array is ${pop_array[*]}"

# echo "pop_array size is ${#pop_array[@]}"

for ((i=0 ; i<${#pop_array[@]} ; i++))
do
    pry_init_1=${pop_array[$i]} # argv[7] prey 1 initial density in nb of individuals
    pry_init_2=$pry_init_1 # argv[8] prey 2 initial density
    prd_init_1=0 # argv[20] predator 1 initial density in nb of individuals

    echo "prey populations initial density are $pry_init_1 and $pry_init_2 and predators $prd_init_1"

    # name the simulation with only the variables of interest and their value
    sim_name="prey1ParamChoice-size$size-res1max$max_res_1-res2max$max_res_2-pry1init$pry_init_1-pry2init$pry_init_2-prdInit$prd_init_1" # argv[1]

    ## Create useful directories ##

    # if the directory exists, delete content? otherwise create it
    while true; do
        if [ -d "$sim_name" ]; then
            read -p "The simulation directory already exists, its content will be deleted, continue?" yn
            case $yn in
                [Yy]* ) find $sim_name -mindepth 1 -delete; break;;
                [Nn]* ) exit;;
                * ) echo "Please answer yes or no.";;
            esac
        else
            # create a directory for the simulation
            mkdir $sim_name
        fi
    done

    # move the executable program there
    cp chapter2ibm.o $sim_name/chapter2ibm.o

    # move to this directory
    cd $sim_name

    ## Start simulation replicates ##

    # launch replicates iteratively
    for ((j=0 ; j<$rep ; j++))
    do
        echo "lauching replicate #$j"

        # seed for random number generator
        rand_seed=$RANDOM # argv[34] set the seed with a random number 

        ./chapter2ibm.o $sim_name $size $res_nb $max_res_1 $max_res_2 $pry_nb $pry_init_1 $pry_init_2 $pry_move_1 $pry_move_2 $pry_cons_1 $pry_cons_2 $pry_surv_1 $pry_surv_2 $pry_offs_1 $pry_offs_2 $pry_repr_1 $pry_repr_2 $prd_nb $prd_init_1 $prd_move_1 $prd_cons_1 $prd_conv_1 $prd_surv_1 $prd_offs_1 $prd_repr_1 $prd_intr_1 $simu_time $burn_time $freq_repr $freq_surv $freq_rslt $freq_snap $rand_seed
        mv $sim_name-ResultsTable.csv rep$j-$sim_name-ResultsTable.csv
        mv $sim_name-SnapshotTable.csv rep$j-$sim_name-SnapshotTable.csv
        # echo "replicate #$j ended"
    done

    echo "All replicates done, moving back to the upper directory"

    # move back to the upper directory for security
    cd ..

done
