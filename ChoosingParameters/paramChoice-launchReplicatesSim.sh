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
#   bash ./paramChoice-launchReplicatesSim.sh


#### Parameter list: assign values ####

## User defined variables ##

# landscape variables
size=4 # argv[2] world's side size
res_nb=2 # argv[3] number of resource types
max_res_1=10 # argv[4] max resource 1 per cell
max_res_2=10 # argv[5] max resource 1 per cell

max_cell=2.5 # max expected number of preys of each kind per cell # if 2.5, between 2 and 3 animals per cell

# prey variables
pry_nb=2 # argv[6] number of prey types
pry_move_1=0.25 # argv[9] prey 1 max movement range in fraction of size
pry_move_2=0.25 # argv[10] prey 2 movement range
pry_offs_1=1 # argv[15] prey 1 max number of offspring
pry_offs_2=1 # argv[16] prey 2 max number of offspring

# predator variables
prd_nb=1 # argv[19] number of predator types
prd_move_1=0.25 # argv[21] predator max movement range in fraction of size
prd_offs_1=1 # argv[23] predator max number of offspring
prd_intr_1=1 # argv[25] predator time of introduction in the model
prd_asym_1=1 # argv[26] asymmetry in prey1 to prey2 conversion rates

# time variables
simu_time=10 # argv[27] simulation time
freq_repr=3 # argv[28] frequency of reproduction trials
freq_surv=3 # argv[29] frequency of survival trials
# burn_time=0 # argv[] time steps to burn before lauching regular survival and reproduction trials

# frequency of assessment
freq_rslt=3 # argv[30] frequency of results snap shot
freq_snap=5 # argv[31] frequency of snapshot measure

# number of replicates
rep=3

## Non user defined variables ##

# pry_cons_1= # argv[11] prey 1 max consumption in resource units
divide=$max_res_1; by=$max_cell; pry_cons_1=`echo "scale=0; ($divide+$by-1)/$by" | bc`; echo "pry_cons_1 = $pry_cons_1" # bash way to ceil a float. by = max expected number of animal per cell HARD CODED NOT IDEAL.
# pry_cons_2= # argv[12] prey 2 max consumption
divide=$max_res_2; by=$max_cell; pry_cons_2=`echo "scale=0; ($divide+$by-1)/$by" | bc`; echo "pry_cons_2 = $pry_cons_2"
# pry_surv_1= # argv[13] prey 1 resource units needed to pass survival trial. 
divide=$freq_surv*$pry_cons_1; by=3; pry_surv_1=`echo "scale=0; ($divide+$by-1)/$by" | bc`; echo "pry_surv_1 = $pry_surv_1" # by = a third of the max number of consecutive fasting days.
# pry_surv_2= # argv[14] prey 2 resource units needed to pass survival trial
divide=$freq_surv*$pry_cons_2; by=3; pry_surv_2=`echo "scale=0; ($divide+$by-1)/$by" | bc`; echo "pry_surv_2 = $pry_surv_2"
pry_repr_1=$pry_surv_1; # argv[17] prey 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.
pry_repr_2=$pry_surv_2; # argv[18] prey 2 resource units needed to pass reproduction trial

prd_surv_1=$(($pry_surv_1*2)); echo "prd_surv_1 = $prd_surv_1" # argv[22] predator1 resource units needed to pass survival trial # defined as a fraction of prey1's
prd_repr_1=$prd_surv_1; echo "prd_repr_1 = $prd_repr_1" # argv[24] predator 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.


#### Simulation loop ####

# max=$size*$size*3*2 # 2 times the max expected number of animal per cell times the number of cells
# min=$size*$size*3/10 # one tenth of the maximum

pop_array=($(seq 10 10 30)) # hard coded for now

# echo "pop_array is ${pop_array[*]}"

# echo "pop_array size is ${#pop_array[@]}"

for ((i=0 ; i<${#pop_array[@]} ; i++))
do
    pry_init_1=${pop_array[$i]} # argv[7] prey 1 initial density in nb of individuals
    # pry_init_1=10 # argv[7] prey 1 initial density in nb of individuals
    pry_init_2=$pry_init_1 # argv[8] prey 2 initial density
    prd_init_1=0 # argv[20] predator 1 initial density in nb of individuals

    echo "prey populations initial density are $pry_init_1 and $pry_init_2 and predators $prd_init_1"

    # name the simulation with only the variables of interest and their value
    sim_name="choosingPreyInitialPop-size$size-res1max$max_res_1-res2max$max_res_2-pry1init$pry_init_1-pry2init$pry_init_2-prdInit$prd_init_1" # argv[1]

    ## Create useful directories ##

    # if the directory exists, delete content? otherwise create it
    while true; do
        if [ -d "$sim_name" ]
        then
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

    # Use make to compile the ibm
    make

    # move the executable program there
    cp chapter2ibm.o $sim_name/chapter2ibm.o

    # move to this directory
    cd $sim_name

    # write a txt file with all the parameter values
    printf "sim_name=$sim_name # argv[1] \n\n" >> paramFile.txt
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
    printf "prd_intr_1=$prd_intr_1 # argv[25] predator 1 time of introduction in the model\n" >> paramFile.txt
    printf "prd_asym_1=$prd_asym_1 # argv[26] predator 1 asymmetry in prey1 to prey2 conversion rates\n\n" >> paramFile.txt
    printf "# time variables\n" >> paramFile.txt
    printf "simu_time=$simu_time # argv[27] simulation time\n" >> paramFile.txt
    printf "freq_repr=$freq_repr # argv[28] frequency of reproduction trials\n" >> paramFile.txt
    printf "freq_surv=$freq_surv # argv[29] frequency of survival trials\n\n" >> paramFile.txt
    printf "# frequency of assessment\n" >> paramFile.txt
    printf "freq_rslt=$freq_rslt # argv[30] frequency of landscape snap shot\n" >> paramFile.txt
    printf "freq_snap=$freq_snap # argv[31] frequency of results measure\n\n" >> paramFile.txt
    printf "freq_snap=$freq_snap # argv[31] frequency of results measure\n\n" >> paramFile.txt
    printf "# number of replicates\n\n" >> paramFile.txt
    printf "rep=$rep\n" >> paramFile.txt

    ## Start simulation replicates ##

    # launch replicates iteratively
    for ((j=0 ; j<$rep ; j++))
    do
        echo "lauching replicate #$j"

        # seed for random number generator
        rand_seed=$RANDOM # argv[32] set the seed with a random number 

        ./chapter2ibm.o $sim_name $size $res_nb $max_res_1 $max_res_2 $pry_nb $pry_init_1 $pry_init_2 $pry_move_1 $pry_move_2 $pry_cons_1 $pry_cons_2 $pry_surv_1 $pry_surv_2 $pry_offs_1 $pry_offs_2 $pry_repr_1 $pry_repr_2 $prd_nb $prd_init_1 $prd_move_1 $prd_surv_1 $prd_offs_1 $prd_repr_1 $prd_intr_1 $prd_asym_1 $simu_time $freq_repr $freq_surv $freq_rslt $freq_snap $rand_seed
        mv $sim_name-ResultsTable.csv rep$j-$sim_name-ResultsTable.csv
        mv $sim_name-SnapshotTable.csv rep$j-$sim_name-SnapshotTable.csv
        # echo "replicate #$j ended"
    done

    echo "All replicates done, moving back to the upper directory"

    # move back to the upper directory before looping
    cd ..

done
