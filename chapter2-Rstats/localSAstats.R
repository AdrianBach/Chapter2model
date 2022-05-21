library(ggplot2)
library(stringr)

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

mergeResults <- function(path, keyword = c("Results", "Snapshot"), pattern) {
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c("folder"), x = content, value = T)
  
  # loop over the folders
  for (i in 1:length(content)) {
    
    # path to folder
    # folder = paste("~/", path, content[i], sep = "")
    folder = paste(path, content[i], sep = "")
    print(paste("in localSA folder ", folder))
    
    # get content
    simFol = list.files(folder)
    
    # select only folders
    simFol <- grep(pattern = c(pattern), x = simFol, value = T)
    
    # loop over the sim folders
    for (j in 1:length(simFol)) {
      
      # get results folder
      resFol <- paste(folder, simFol[j], sep = "/")
      print(paste("in sim folder ", resFol))
      
      # get files
      results <- list.files(resFol)
      
      # only csv files
      results <- grep(pattern = c(".csv"), x = results, value = T)
      
      # only the results files
      results <- grep(pattern = c(keyword), x = results, value = T)
      
      # loop over the results files
      for (m in 1:length(results)) {
        
        # print(paste("reading file ", results[m]))
        # read results
        res <- read.csv(paste(resFol, results[m], sep = "/"))
          
        # if first round, create merge table
        if (m == 1) {
          
          # create an empty table
          headers <- colnames(res)
          headers <- c("replicate", colnames(res), "prey1growthRate", "prey2growthRate", "predatorGrowthRate")
          tab <- data.frame(matrix(ncol = length(headers), nrow = 0))
        }
        
        # bind the replicate number to the table
        replicate <- rep(m, times = dim(res)[1])
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
      newFolderDir <- paste(resFol, "/stats-", simFol[j], sep = "")
      dir.create(path = newFolderDir)
      write.csv(tab, file = paste(newFolderDir, "/merged", keyword, "-", simFol[j], sep = ""), row.names = FALSE)
    
    } # end loop over sim folders
    
  } # end loop over folders
  
} # end of function

statsResults <- function(path, keyword = c("Results", "Snapshot"), pattern) {
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c("folder"), x = content, value = T)
  
  # loop over the local SA folders
  for (i in 1:length(content)) {
    
    # path to folder
    folder = paste(path, content[i], sep = "")
    print(paste("in localSA folder ", folder))
    
    # create a global stats folder
    dir.create(path = paste(folder, "allStatsAndPlots", sep = "/"))
    
    # create a Results folder
    dir.create(path = paste(folder, "/allStatsAndPlots/", keyword, "Files", sep = ""))
    
    # get content
    simFol = list.files(folder)
    
    # select only folders
    simFol <- grep(pattern = c(pattern), x = simFol, value = T)
    
    # loop over the sim folders
    for (j in 1:length(simFol)) {
      
      print(paste("in sim folder ", simFol[j]))
      
      # # rename the file
      # cux <- list.files(paste(path, content[i], simFol[j], sep = "/"))
      # cux <- grep(pattern = c("stats"), x = cux, value = T)
      # file.rename(paste(path, content[i], simFol[j], cux, sep = "/"), paste(path, content[i], "/", simFol[j], "/stats-", simFol[j], sep = ""))
      
      # path to folder
      statsFol = paste(path, content[i], "/", simFol[j], "/stats-", simFol[j], sep = "")
      
      # file.rename(paste(statsFol, results, sep = "/"), paste(statsFol, "/merged", keyword, "-", simFol[j], ".csv", sep = ""))
      
      # get files
      results <- list.files(statsFol)
      
      # only csv files
      results <- grep(pattern = c("merged"), x = results, value = T)
      
      # only the results files
      results <- grep(pattern = c(keyword), x = results, value = T)
      
      # get name
      name <- paste(statsFol, "/", results, sep = "")
      
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
      for (k in 1:length(ts)) {
        
        # subset per ts number
        sub <- subset(stats, stats$timeStep == as.numeric(ts[k]))
        
        # initiate new line with time step
        newLine <- c(as.numeric(ts[k]), dim(sub)[1])
        
        # apply stats to each column of measures
        for (k in 3:(dim(sub)[2])) { # first line is rep, second is ts
          newLine <- c(newLine, mean(sub[,k]), boot_sd_ci(sub[,k])[2], boot_sd_ci(sub[,k])[3])
        }
        
        # rbind the line to tab
        tab <- rbind(tab, newLine)
      } # end loop over time steps
      
      # column names
      colnames(tab) <- newHeaders
      
      # write in the corresponding folder
      write.csv(tab, file = paste(statsFol, "/stats", keyword, "-", simFol[j], ".csv", sep = ""), row.names = FALSE)
      
      # and copy in the global stats folder
      file.copy(from = paste(statsFol, "/stats", keyword, "-", simFol[j], ".csv", sep = ""), to = paste(folder, "/allStatsAndPlots/", keyword, "Files", sep = ""))
      
      # plot and save figure
      data <- tab
      
      x = data$timeStep
      y1 = data$prey1PopulationSizeMean
      y2 = data$prey2PopulationSizeMean
      y3 = data$predator1PopulationSizeMean
      y1min = data$prey1PopulationSizeICinf
      y2min = data$prey2PopulationSizeICinf
      y3min = data$predator1PopulationSizeICinf
      y1max = data$prey1PopulationSizeICsup
      y2max = data$prey2PopulationSizeICsup
      y3max = data$predator1PopulationSizeICsup
      y1c = "red"
      y2c = "blue"
      y3c = "orange"
      tIntro = 210
      
      fig <- ggplot(data, aes(x)) + 
        geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
        geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
        geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
        geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
        geom_line(aes(y = y1), color = y1c) +
        geom_line(aes(y = y2), color = y2c) +
        geom_line(aes(y = y3), color = y3c) +
        # geom_point(aes(y = y1), size = 2.5, shape = 21, fill = "white", color = y1c) +
        # geom_point(aes(y = y2), size = 2.5, shape = 22, fill = "white", color = y2c) +
        # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c) +
        labs(x = "Time steps", y = "Population size") +
        scale_colour_manual(name='Populations',
                            breaks=c('Prey 1', 'Prey 2', 'Predator'),
                            values=c(y1c, y2c, y3c))
      
      # save plot in this folder
      ggsave(filename = paste("stats", keyword, "-", simFol[j], ".pdf", sep = ""), path = statsFol, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
      
      # and copy in the global stats folder
      file.copy(from = paste(statsFol, "/stats", keyword, "-", simFol[j], ".pdf", sep = ""), to = paste(folder, "/allStatsAndPlots/", keyword, "Files", sep = ""))
      
    } # end loop over sim folders
    
  } # end loop over local SA folders
  
} # end of function

localSAres <- function(path, keyword = c("Results", "Snapshot"), pattern) {
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c("folder"), x = content, value = T)
  
  # loop over the folders
  for (i in 1:length(content)) {
    
    # path to folder
    # folder = paste("~/", path, content[i], sep = "")
    folder = paste(path, content[i], sep = "")
    print(paste("in localSA folder ", folder))
    
    # get content
    simFol = list.files(folder)
    
    # select only folders
    simFol <- grep(pattern = c(pattern), x = simFol, value = T)
    
    # loop over the sim folders
    for (j in 1:length(simFol)) {
      
      # get results folder
      resFol <- paste(folder, simFol[j], sep = "/")
      print(paste("in sim folder ", resFol))
      
      # get files
      results <- list.files(resFol)
      
      # only csv files
      results <- grep(pattern = c(".csv"), x = results, value = T)
      
      # only the results files
      results <- grep(pattern = c(keyword), x = results, value = T)
      
      # loop over the results files
      for (m in 1:length(results)) {
        
        # print(paste("reading file ", results[m]))
        # read results
        res <- read.csv(paste(resFol, results[m], sep = "/"))
        
        # count extinctions
        
      } # end loop over files (replicates)
      
      # name columns
      colnames(tab) <- headers
      
      # create stats folder and save table
      newFolderDir <- paste(resFol, "/stats-", simFol[j], sep = "")
      dir.create(path = newFolderDir)
      write.csv(tab, file = paste(newFolderDir, "/merged", keyword, "-", simFol[j], sep = ""), row.names = FALSE)
      
    } # end loop over sim folders
    
  } # end loop over folders
  
}

Path = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/"
Pattern = "localSA-pred"
Keyword = "Results"

mergeResults(path = Path, keyword = Keyword, pattern = Pattern)
statsResults(path = Path, keyword = Keyword, pattern = Pattern)

