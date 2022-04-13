boot_sd_ci <- function(x, confidence = 95, itr = 1000) {
  
  # init iterations and sample
  i <- itr
  bt_avg <- NULL
  
  # loop over iterations
  while (i > 0) {
    # sample randomly from x
    spl <- sample(x, length(x), replace = TRUE)
    
    # store mean
    bt_avg <- c(bt_avg, sum(spl)/length(spl))
    
    # decrement i
    i <- i-1
  }
  
  # mean over the bootstrapped samples
  bt_est <- sum(bt_avg)/itr
  
  # compute standard deviation
  bt_sd <- sqrt((1/(length(x)-1)) * sum((bt_avg-bt_est)^2))
  
  # compute confidence interval
  # sorting bt_avg numerically
  st_avg <- sort(bt_avg)
  
  # get the first value after the 2.5 first centiles
  bt_95ci_inf <- st_avg[floor(0.5*(1-0.01*confidence)*itr)+1]
  
  # get the last value before the 2.5 last centiles
  bt_95ci_sup <- st_avg[floor((0.01*confidence+0.5*(1-0.01*confidence))*itr)-1]
  
  res <- c(bt_sd, bt_95ci_inf, bt_95ci_sup)
  return(res) 
  
}

mergeResults <- function(path, keyword = c("Results", "Snapshot"), Pattern) {
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = Pattern, x = content, value = T)
  
  # loop over the folders
  for (i in 1:length(content)) {
    
    # path to folder
    # folder = paste("~/", path, content[i], sep = "")
    folder = paste(path, content[i], sep = "/")
    
    # get content
    sim = list.files(folder)
    
    # select only csv files
    results <- grep(pattern = c(".csv"), x = sim, value = T)
    
    # select only results files
    results <- grep(pattern = c(keyword), x = results, value = T)
    
    # loop over the files
    for (j in 1:length(results)) {
      
      # read one file to get column number and headers
      res <- read.csv(paste(folder, results[j], sep = "/"))
      
      # if first round, create merge table
      if (j == 1) {
              
        # create an empty table
        headers <- colnames(res)
        headers <- c("replicate", colnames(res), "prey1growthRate", "prey2growthRate", "predatorGrowthRate")
        tab <- data.frame(matrix(ncol = length(headers), nrow = 0))
      }
      
      # bind the replicate number to the table
      replicate <- rep(j, times = dim(res)[1])
      res <- cbind(replicate, res)
      
      # Compute growth rate
      prey1growthRate <- 0
      prey2growthRate <- 0
      predatorGrowthRate <- 0
      
      # loop over the lines
      for (k in 2:dim(res)[1]) {
        prey1growthRate     <- c(prey1growthRate, ifelse(res$prey1PopulationSize[k-1] != 0, 100*(res$prey1PopulationSize[k]-res$prey1PopulationSize[k-1])/res$prey1PopulationSize[k-1], 0))
        prey2growthRate     <- c(prey2growthRate, ifelse(res$prey2PopulationSize[k-1] != 0, 100*(res$prey2PopulationSize[k]-res$prey2PopulationSize[k-1])/res$prey2PopulationSize[k-1], 0))
        predatorGrowthRate  <- c(predatorGrowthRate, ifelse(res$predator1PopulationSize[k-1] != 0, 100*(res$predator1PopulationSize[k]-res$predator1PopulationSize[k-1])/res$predator1PopulationSize[k-1], 0))
      } # end loop over lines
      
      res <- cbind(res, prey1growthRate, prey2growthRate, predatorGrowthRate)
      
      tab <- rbind(tab, res)
    } # end loop over files (replicates)
    
    # name columns
    colnames(tab) <- headers
    
    # create stats folder and save table
    newFolderDir <- paste(folder, "/stats-", content[i], sep = "")
    dir.create(path = newFolderDir)
    write.csv(tab, file = paste(newFolderDir, "/merged", keyword, "-", content[i], '.csv', sep = ""), row.names = FALSE)
  } # end loop over folders
  
} # end of function

statsResults <- function(path, keyword = c("Results", "Snapshot"), Pattern) {
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c(Pattern), x = content, value = T)

  # loop over the folders
  for (i in 1:length(content)) {
    
    # path to folder
    # folder = paste("~/", path, content[i], "/stats-", content[i], sep = "")
    folder = paste(path, "/", content[i], "/stats-", content[i], sep = "")
    
    # get content
    sim = list.files(folder)
    
    # select only csv files
    results <- grep(pattern = c("merged"), x = sim, value = T)
    
    # select only results files
    results <- grep(pattern = c(keyword), x = results, value = T)
    
    # get name
    name <- paste(folder, "/", results, sep = "")
    
    # read table
    stats <- read.csv(name)
    
    ## create table
    
    # new headers
    headers <- colnames(stats)
    newHeaders <- c(headers[2], "repNb")
    for (k in 3:length(headers)) {
      newHeaders <- c(newHeaders, paste(headers[k], "Mean", sep = ""), paste(headers[k], "ICinf", sep = ""), paste(headers[k], "ICsup", sep = ""))
    }
    
    # create an empty table 
    tab <- data.frame(matrix(ncol = length(newHeaders), nrow = 0))
      
    # for loop making subset for each time step
    ts <- levels(as.factor(stats$timeStep))
    for (j in 1:length(ts)) {
      
      # subset per ts number
      sub <- subset(stats, stats$timeStep == as.numeric(ts[j]))
      
      # initiate new line with time step
      newLine <- c(as.numeric(ts[j]), dim(sub)[1])
      
      # apply stats to each column of measures
      for (k in 3:(dim(sub)[2])) { # first line is rep, second is ts
        newLine <- c(newLine, mean(sub[,k]), boot_sd_ci(sub[,k])[2], boot_sd_ci(sub[,k])[3])
      }
      
      # rbind the line to tab
      tab <- rbind(tab, newLine)
    }
    
    # column names
    colnames(tab) <- newHeaders
    
    # write in the corresponding folder
    write.csv(tab, file = paste(folder, "/stats", keyword, "-", content[i], ".csv", sep = ""), row.names = FALSE)
    
  } # end of loop over folders
} # end of function

getPosStab <- function(array, threshold, period) {
  
  posStab <- 1
  
  # initialise count
  count <- 0
  
  # parse array
  for (i in 1:length(array)) {
    if (abs(array[i]) < threshold && abs(array[i]) > 0) {
      count = count+1
    } else {
      count = 0
    } 
    
    if (count >= period) {
      posStab = i
      break
    }
  }
  
  return(posStab)
}

posStabNoPred <- function(path, keyword = c("Results", "Snapshot"), Pattern, stabThreshold, stabPeriod) {
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c(Pattern), x = content, value = T)
    
  ## create table
  # new headers
  headers <- c("preyInit", "repNb", "threshold", "period", "pry1posStab", "pry2posStab", "prdPosStab", "prey1popEq", "prey2popEq", "predatorPopEq")
  newHeaders <- headers[1:(length(headers)-3)]
  for (k in (length(headers)-2):length(headers)) {
    newHeaders <- c(newHeaders, paste(headers[k], "Mean", sep = ""), paste(headers[k], "ICinf", sep = ""), paste(headers[k], "ICsup", sep = ""))
  }
  
  # create an empty table 
  tab <- data.frame(matrix(ncol = length(newHeaders), nrow = 0))
  
  # loop over the folders
  for (i in 1:length(content)) {
    
    # path to folder
    # folder = paste("~/", path, content[i], "/stats-", content[i], sep = "")
    folder = paste(path, "/", content[i], "/stats-", content[i], sep = "")
    
    # get content
    results = list.files(folder)
    
    # select only csv files
    results <- grep(pattern = c("stats"), x = results, value = T)
    
    # select only results files
    results <- grep(pattern = c(keyword), x = results, value = T)
    
    # get name
    name <- paste(folder, results, sep = "/")
    
    # read table
    stats <- read.csv(name)
    
    # find growthRate stabilization time step
    pry1tstab <- getPosStab(stats$prey1growthRateMean, stabThreshold, stabPeriod)
    pry2tstab <- getPosStab(stats$prey2growthRateMean, stabThreshold, stabPeriod)
    prdTstab <- getPosStab(stats$predatorGrowthRateMean, stabThreshold, stabPeriod)
    
    newLine <- c(stats$prey1PopulationSizeMean[1], stats$repNb[1], stabThreshold, stabPeriod, stats$timeStep[pry1tstab], stats$timeStep[pry2tstab], stats$timeStep[prdTstab], 
                 stats$prey1PopulationSizeMean[pry1tstab], stats$prey1PopulationSizeICinf[pry1tstab], stats$prey1PopulationSizeICsup[pry1tstab], 
                 stats$prey2PopulationSizeMean[pry2tstab], stats$prey2PopulationSizeICinf[pry2tstab], stats$prey2PopulationSizeICsup[pry2tstab], 
                 stats$predator1PopulationSizeMean[prdTstab], stats$predator1PopulationSizeICinf[prdTstab], stats$predator1PopulationSizeICsup[prdTstab])
      
    # rbind the line to tab
    tab <- rbind(tab, newLine)
  } # end of loop over folders
    
  # column names
  colnames(tab) <- newHeaders
  
  # write in the corresponding folder
  write.csv(tab, file = paste(folder, "equilibriumNoPred.csv", sep = "/"), row.names = FALSE)

} # end of function

# path = paste(getwd(), "Duthie_Sims", sep="/")
path = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/"
# path = getwd()
keyword = "Results"
Pattern = "test-size"

mergeResults(path, keyword, Pattern)
statsResults(path, keyword, Pattern)

stabPeriod = 5
stabThreshold = 10

posStabNoPred(path, keyword, Pattern, stabThreshold, stabPeriod)
