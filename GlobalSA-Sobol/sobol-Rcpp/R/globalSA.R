# Adrian BACH
# script to run Sobol method
# guidelines: https://www.jasss.org/17/3/11.html §3.38 - §3.46

library(Rcpp)

path = getwd()

######## list of parameters ########

#### User defined variables ####

{
# simulation name
sim_name = "globalSA-sobol-sim"
  
# landscape variables
size=25     # argv[2] world's side size
res_nb=2     # argv[3] number of resource types
max_res_1=50 # argv[4] max resource 1 per cell
max_res_2=50 # argv[5] max resource 1 per cell

max_cell=3.3 # max expected number of preys of each kind per cell # if 2.5, between 2 and 3 animals per cell

# prey variables
pry_nb=2        # argv[6] number of prey types
pry_init_1=250  # argv[7] prey 1 initial density in nb of individuals
pry_init_2=pry_init_1 # argv[8] prey 2 initial density
pry_move_1=0.1  # argv[9] prey 1 max movement range in fraction of size
pry_move_2=0.1  # argv[10] prey 2 movement range
pry_offs_1=1    # argv[15] prey 1 max number of offspring
pry_offs_2=1    # argv[16] prey 2 max number of offspring

# predator variables
prd_nb=1        # argv[19] number of predator types
prd_init_1=25   # argv[20] predator 1 initial density in nb of individuals
prd_move_1=0.1  # argv[21] predator 1 max movement range in fraction of size
prd_offs_1=1    # argv[24] predator 1 max number of offspring
prd_intr_1=201  # argv[26] predator 1 time of introduction in the model
prd_ctch_pry1_1=0.1  # argv[27] predator 1 prey1 catch probability
prd_ctch_pry2_1=0.1  # argv[28] predator 1 prey2 catch probability
prd_oprt_1=0    # argv[31] is predator oportunistic? (0 or 1)
prd_spcf_1=0    # argv[32] is predator specific? (0 or 1)

# time variables
simu_time=1000   # argv[33] simulation time
freq_surv=10     # argv[35] frequency of survival trials
freq_repr=freq_surv     # argv[34] frequency of reproduction trials
freq_rfll=freq_surv     # argv[36] frequency of landscape resources refill

# frequency of assessment
freq_rslt=freq_surv    # argv[37] frequency of landscape results shot
freq_snap=100   # argv[38] frequency of snapshot measure

#### Non user defined variables ####

# preys variables
pry_cons_1 = ceiling(max_res_1/max_cell) # argv[11] prey 1 max consumption in resource units
pry_cons_2 = ceiling(max_res_2/max_cell) # argv[12] prey 2 max consumption
pry_surv_1 = ceiling(freq_surv*pry_cons_1/3) # argv[13] prey 1 resource units needed to pass survival trial. 
pry_surv_2 = ceiling(freq_surv*pry_cons_2/3) # argv[14] prey 2 resource units needed to pass survival trial
pry_repr_1=pry_surv_1; # argv[17] prey 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.
pry_repr_2=pry_surv_2; # argv[18] prey 2 resource units needed to pass reproduction trial

# predator variables
prd_cons_1=3*pry_cons_1 # arg[22]
prd_surv_1=ceiling(prd_cons_1*freq_surv/3) # arg[23]
prd_repr_1=prd_surv_1 # argv[25] predator 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.
prd_cvrt_pry1_1=floor(freq_surv*prd_cons_1/3) # argv[29] predator 1 prey1 resources/catch
ratio=1
prd_cvrt_pry2_1=ratio*prd_cvrt_pry1_1  # argv[30] predator 1 prey1 resources/catch
}

# store in a vector to "replace" argv ?
params <- c(size,             # argv[2] world's side size
            res_nb,           # argv[3] number of resource types
            max_res_1,        # argv[4] max resource 1 per cell
            max_res_2,        # argv[5] max resource 1 per cell
            pry_nb,           # argv[6] number of prey types
            pry_init_1,       # argv[7] prey 1 initial density in nb of individuals
            pry_init_2,       # argv[8] prey 2 initial density
            pry_move_1,       # argv[9] prey 1 max movement range in fraction of size
            pry_move_2,       # argv[10] prey 2 movement range
            pry_cons_1,       # argv[11] prey 1 max consumption in resource units
            pry_cons_2,       # argv[12] prey 2 max consumption
            pry_surv_1,       # argv[13] prey 1 resource units needed to pass survival trial. 
            pry_surv_2,       # argv[14] prey 2 resource units needed to pass survival trial
            pry_offs_1,       # argv[15] prey 1 max number of offspring
            pry_offs_2,       # argv[16] prey 2 max number of offspring
            pry_repr_1,       # argv[17] prey 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.
            pry_repr_2,       # argv[18] prey 2 resource units needed to pass reproduction trial
            prd_nb,           # argv[19] number of predator types
            prd_init_1,       # argv[20] predator 1 initial density in nb of individuals
            prd_move_1,       # argv[21] predator 1 max movement range in fraction of size
            prd_cons_1,       # argv[22]
            prd_surv_1,       # argv[23]
            prd_offs_1,       # argv[24] predator 1 max number of offspring
            prd_repr_1,       # argv[25] predator 1 resource units needed to pass reproduction trial. Defined as a proportion of what is needed to pass survival trial.
            prd_intr_1,       # argv[26] predator 1 time of introduction in the model
            prd_ctch_pry1_1,  # argv[27] predator 1 prey1 catch probability
            prd_ctch_pry2_1,  # argv[28] predator 1 prey2 catch probability
            prd_cvrt_pry1_1,  # argv[29] predator 1 prey1 resources/catch
            prd_cvrt_pry2_1,  # argv[30] predator 1 prey1 resources/catch
            prd_oprt_1,       # argv[31] is predator oportunistic? (0 or 1)
            prd_spcf_1,       # argv[32] is predator specific? (0 or 1)
            simu_time,        # argv[33] simulation time
            freq_surv,        # argv[35] frequency of survival trials
            freq_repr,        # argv[34] frequency of reproduction trials
            freq_rfll,        # argv[36] frequency of landscape resources refill
            freq_rslt,        # argv[37] frequency of landscape results shot
            freq_snap)        # argv[38] frequency of snapshot measure

#### function inputing paramaters outputting measure ####

sim <- function(params) {
  
  ## simulation name
  # sim_name=paste("globalSA-predSpcf", params[32], "-predOpnt", params[31], "-pryConvRateRatio",params[30]/prd_cvrt_pry1_1[29], "-pryCtchProbRatio1", params[28]/params[27], "-pryOfspRatio1", params[16]/params[15],"-pryMaxConsRatio1", params[12]/params[11], sep = "") # argv[1]
  # or a simpler name defined earlier to be written over everytime, might be better for data load
  
  ## add as params[1]
  params <- c(sim_name, params)
  
  ## run model with this params vector
  sourceCpp(paste(path, '/src/chapter2model-v0.5.0-Rrun.cpp', sep = ""))
  tail(main(params)) # will params replace the argv vector???
  
  ## Try to get one value to feed criterias and see if it works
  
  # read in the results file
  # the file name should be: paste(params[1], "-Results.csv", sep = "")
  
  # get the final value of prey 1 density for example
  
  # criterias =
  
  ## compute all criterias
  # inspiration from localSA-stats
  
  ## c them on one vector
  # criterias <- c()
    
  return(criterias)
}

######## Generate two random LHC sets ########

######## Create a Sobol instance, whatever that means ########

######## Add results to the Sobol instance, whatever that means ########