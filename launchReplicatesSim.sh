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
size=4 # argv[2] world's side size
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
pry_offs_1=2 # argv[15] prey 1 max number of offspring
pry_offs_2=2 # argv[16] prey 2 max number of offspring
pry_repr_1=5 # argv[17] prey 1 resource units needed to pass reproduction trial
pry_repr_2=5 # argv[18] prey 2 resource units needed to pass reproduction trial

# predator variables
prd_nb=1 # argv[19] number of predator types
prd_init_1=5 # argv[20] predator 1 initial density in nb of individuals
prd_move_1=0.3 # argv[21] predator 1 max movement range in fraction of size
prd_cons_1=1 # argv[22] predator 1 max catches per time step in prey unit
prd_conv_1=5 # argv[23] predator 1 conversion rate of 1 catch into resources
prd_surv_1=3 # argv[24] predator 1 resource units needed to pass survival trial
prd_offs_1=2 # argv[25] predator 1 max number of offspring
prd_repr_1=5 # argv[26] predator 1 resource units needed to pass reproduction trial
prd_intr_1=4 # argv[27] predator 1 time of introduction in the model

# time variables
simu_time=20 # argv[28] simulation time
burn_time=0 # argv[29] time steps to burn before lauching regular survival and reproduction trials
freq_repr=2 # argv[30] frequency of reproduction trials
freq_surv=2 # argv[31] frequency of survival trials

# frequency of assessment
freq_rslt=1 # argv[32] frequency of landscape snap shot
freq_snap=5 # argv[33] frequency of results measure

# name the simulation with only the variables of interest and their value
sim_name="control-size$size-res1max$max_res_1-res2max$max_res_2-pry1init$pry_init_1-pry2init$pry_init_2-prdInit$prd_init_1" # argv[1]

# number of replicates
rep=3

#### Create useful directories ####

# if the directory exists, delete content? otherwise create it
while true; do
    if [ -d $sim_name ]; then
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
mv chapter2ibm.o $sim_name/chapter2ibm.o

# move to this directory
cd $sim_name

# write a txt file with all the parameter values
printf "sim_name=$sim_name # argv[1] \n" >> paramFile.txt
printf "# landscape variables\n" >> paramFile.txt
printf "size=$size # argv[2] world's side size\n" >> paramFile.txt
printf "res_nb=$res_nb # argv[3] number of resource types\n" >> paramFile.txt
printf "max_res_1=$max_res_1 # argv[4] max resource 1 per cell\n" >> paramFile.txt
printf "max_res_2=$max_res_2 # argv[5] max resource 1 per cell\n\n" >> paramFile.txt
printf "# prey variables\n" >> paramFile.txt
printf "pry_nb=$pry_nb # argv[6] number of prey types\n" >> paramFile.txt
printf "pry_init_1=$pry_init_1 # argv[7] prey 1 initial density in nb of individuals\n" >> paramFile.txt
printf "pry_init_2=$pry_init_2 # argv[8] prey 2 initial density\n" >> paramFile.txt
printf "pry_move_1=$pry_move_1 # argv[9] prey 1 max movement range in fraction of size\n" >> paramFile.txt
printf "pry_move_2=$pry_move_2 # argv[10] prey 2 movement range\n" >> paramFile.txt
printf "pry_cons_1=$pry_cons_1 # argv[11] prey 1 max consumption in resource units\n" >> paramFile.txt
printf "pry_cons_2=$pry_cons_2 # argv[12] prey 2 max consumption\n" >> paramFile.txt
printf "pry_surv_1=$pry_surv_1 # argv[13] prey 1 resource units needed to pass survival trial\n" >> paramFile.txt
printf "pry_surv_2=$pry_surv_2 # argv[14] prey 2 resource units needed to pass survival trial\n" >> paramFile.txt
printf "pry_offs_1=$pry_offs_1 # argv[15] prey 1 max number of offspring\n" >> paramFile.txt
printf "pry_offs_2=$pry_offs_2 # argv[16] prey 2 max number of offspring\n" >> paramFile.txt
printf "pry_repr_1=$pry_repr_1 # argv[17] prey 1 resource units needed to pass reproduction trial\n" >> paramFile.txt
printf "pry_repr_2=$pry_repr_2 # argv[18] prey 2 resource units needed to pass reproduction trial\n\n" >> paramFile.txt
printf "# predator variables\n" >> paramFile.txt
printf "prd_nb=$prd_nb # argv[19] number of predator types\n" >> paramFile.txt
printf "prd_init_1=$prd_init_1 # argv[20] predator 1 initial density in nb of individuals\n" >> paramFile.txt
printf "prd_move_1=$prd_move_1 # argv[21] predator 1 max movement range in fraction of size\n" >> paramFile.txt
printf "prd_surv_1=$prd_surv_1 # argv[22] predator 1 resource units needed to pass survival trial\n" >> paramFile.txt
printf "prd_offs_1=$prd_offs_1 # argv[23] predator 1 max number of offspring\n" >> paramFile.txt
printf "prd_repr_1=$prd_repr_1 # argv[24] predator 1 resource units needed to pass reproduction trial\n" >> paramFile.txt
printf "prd_intr_1=$prd_intr_1 # argv[25] predator 1 time of introduction in the model\n\n" >> paramFile.txt
printf "# time variables\n" >> paramFile.txt
printf "simu_time=$simu_time # argv[26] simulation time\n" >> paramFile.txt
printf "burn_time=$burn_time # argv[27] time steps to burn before lauching regular survival and reproduction trials\n" >> paramFile.txt
printf "freq_repr=$freq_repr # argv[28] frequency of reproduction trials\n" >> paramFile.txt
printf "freq_surv=$freq_surv # argv[29] frequency of survival trials\n\n" >> paramFile.txt
printf "# frequency of assessment\n" >> paramFile.txt
printf "freq_rslt=$freq_rslt # argv[30] frequency of landscape snap shot\n" >> paramFile.txt
printf "freq_snap=$freq_snap # argv[31] frequency of results measure\n\n" >> paramFile.txt
printf "# name the simulation with only the variables of interest and their value\n" >> paramFile.txt
printf "# number of replicates\n" >> paramFile.txt
printf "rep=$rep\n" >> paramFile.txt

#### Start simulation resplicates ####

# launch replicates iteratively
# for ((i=0 ; i<$rep ; i++))
# do
#     echo "lauching replicate #$i"

#     # seed for random number generator
#     rand_seed=$RANDOM # argv[34] set the seed with a random number 

#     ./chapter2ibm.o $sim_name $size $res_nb $max_res_1 $max_res_2 $pry_nb $pry_init_1 $pry_init_2 $pry_move_1 $pry_move_2 $pry_cons_1 $pry_cons_2 $pry_surv_1 $pry_surv_2 $pry_offs_1 $pry_offs_2 $pry_repr_1 $pry_repr_2 $prd_nb $prd_init_1 $prd_move_1 $prd_cons_1 $prd_conv_1 $prd_surv_1 $prd_offs_1 $prd_repr_1 $prd_intr_1 $simu_time $burn_time $freq_repr $freq_surv $freq_rslt $freq_snap $rand_seed
#     mv $sim_name-ResultsTable.csv rep$i-$sim_name-ResultsTable.csv
#     mv $sim_name-SnapshotTable.csv rep$i-$sim_name-SnapshotTable.csv
#     # echo "replicate #$i ended"
# done

echo "All replicates done, moving back to the upper directory"

# move back to the upper directory for security
cd ..
