#!/bin/bash

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

# landscape variables
size=3 # argv[2] world's side size
res_nb=2 # argv[3] number of resource types
max_res_1=10 # argv[4] max resource 1 per cell
max_res_2=10 # argv[5] max resource 1 per cell

# prey variables
pry_nb=2 # argv[6] number of prey types
pry_init_1=10 # argv[7] prey 1 initial density in nb of individuals
pry_init_2=10 # argv[8] prey 2 initial density
pry_move_1=0.3 # argv[9] prey 1 max movement range in fraction of size
pry_move_2=0.3 # argv[10] prey 2 movement range
pry_cons_1=5 # argv[11] prey 1 max consumption in resource units
pry_cons_2=5 # argv[12] prey 2 max consumption
pry_surv_1=3 # argv[13] prey 1 resource units needed to pass survival trial
pry_surv_2=3 # argv[14] prey 2 resource units needed to pass survival trial
pry_offs_1=3 # argv[15] prey 1 max number of offspring
pry_offs_2=3 # argv[16] prey 2 max number of offspring
pry_repr_1=6 # argv[17] prey 1 resource units needed to pass reproduction trial
pry_repr_2=6 # argv[18] prey 2 resource units needed to pass reproduction trial

# predator variables
prd_nb=1 # argv[19] number of predator types
prd_init_1=5 # argv[20] predator 1 initial density in nb of individuals
prd_move_1=0.3 # argv[21] predator 1 max movement range in fraction of size
prd_cons_1=1 # argv[22] predator 1 max catches per time step in prey unit
prd_conv_1=5 # argv[23] predator 1 conversion rate of 1 catch into resources
prd_surv_1=5 # argv[24] predator 1 resource units needed to pass survival trial
prd_offs_1=2 # argv[25] predator 1 max number of offspring
prd_repr_1=6 # argv[26] predator 1 resource units needed to pass reproduction trial
prd_intr_1=1 # argv[27] predator 1 time of introduction in the model

# time variables
simu_time=4 # argv[28] simulation time
burn_time=1 # argv[29] time steps to burn before lauching regular survival and reproduction trials
freq_repr=1 # argv[30] frequency of reproduction trials
freq_surv=1 # argv[31] frequency of survival trials

# frequency of assessment
freq_snap=1 # argv[32] frequency of landscape snap shot
freq_rslt=1 # argv[33] frequency of results measure

# name the simulation with only the variables of interest and their value
sim_name="control-size$size-res1max$max_res_1-res2max$max_res_2-pry1init$pry_init_1-pry2init$pry_init_2-prdInit$prd_init_1" # argv[1]

# number of replicates
rep=3

#### Create useful directories ####

# create a directory for the simulation
mkdir $sim_name

# copy the executable program there
cp chapter2ibm.o $sim_name/chapter2ibm.o

# move to this directory
cd $sim_name

#### Start simulation resplicates ####

# launch replicates iteratively
for ((i=0 ; i<$rep ; i++))
do
    echo "lauching replicate #$i"
    ./chapter2ibm.o $sim_name $size $res_nb $max_res_1 $max_res_2 $pry_nb $pry_init_1 $pry_init_2 $pry_move_1 $pry_move_2 $pry_cons_1 $pry_cons_2 $pry_surv_1 $pry_surv_2 $pry_offs_1 $pry_offs_2 $pry_repr_1 $pry_repr_2 $prd_nb $prd_init_1 $prd_move_1 $prd_cons_1 $prd_conv_1 $prd_surv_1 $prd_offs_1 $prd_repr_1 $prd_intr_1 $simu_time $burn_time $freq_repr $freq_surv
    mv $sim_name-ResultsTable.csv rep$i-$sim_name-ResultsTable.csv
    mv $sim_name-SnapshotTable.csv rep$i-$sim_name-SnapshotTable.csv
    # echo "replicate #$i ended"
done

echo "All replicates done, moving back to the upper directory"

# move back to the upper directory for security
cd ..
