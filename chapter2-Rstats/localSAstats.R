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
      write.csv(tab, file = paste(newFolderDir, "/merged", keyword, "-", simFol[j], ".csv", sep = ""), row.names = FALSE)
    
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

localSAresults <- function(path, keyword = c("Results", "Snapshot"), pattern) {
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c("folder"), x = content, value = T)
  
  # create table 
  headers <- c("predSpecific", "predOportunistic", "convRateRatio", "catchProbaRatio", "preyOffspRatio", "maxConsRatio", "replicatesNb",
               "prey1extFreq", "prey2extFreq", "pred1extFreq",
               "prey1densBeforeMean", "prey1densBeforeMax", "prey1densBeforeMin",
               "prey2densBeforeMean", "prey2densBeforeMax", "prey2densBeforeMin",
               "prey1growthBeforeMean", "prey1growthBeforeMax", "prey1growthBeforeMin",
               "prey2growthBeforeMean", "prey2growthBeforeMax", "prey2growthBeforeMin",
               "prey1densAfterMean", "prey1densAfterMax", "prey1densAfterMin",
               "prey2densAfterMean", "prey2densAfterMax", "prey2densAfterMin",
               "predatorDensMean", "predatorDensMax", "predatorDensMin",
               "prey1growthAfterMean", "prey1growthAfterMax", "prey1growthAfterMin",
               "prey2growthAfterMean", "prey2growthAfterMax", "prey2growthAfterMin",
               "predatorGrowthMean", "predatorGrowthMax", "predatorGrowthMin",
               "prey1catchesMean", "prey1catchesMax", "prey1catchesMin",
               "prey2catchesMean", "prey2catchesMax", "prey2catchesMin")
  
  # set before after intervals
  tIntro = 210
  tEnd = 1000
  before <- c(tIntro-100, tIntro)
  after <- c(tEnd-250, tEnd)
  
  # loop over the folders
  for (i in 1:length(content)) {
    
    # path to folder
    # folder = paste("~/", path, content[i], sep = "")
    folder = paste(path, content[i], sep = "")
    print(paste("in localSA folder ", i, ": ", folder))# get the value of XparamName
    
    tab <- data.frame(matrix(ncol = length(headers), nrow = 0))
    
    # # name columns
    # colnames(tab) <- headers
    
    # get content
    simFol = list.files(folder)
    
    # select only folders
    simFol <- grep(pattern = c(pattern), x = simFol, value = T)
    
    # loop over the sim folders
    for (j in 1:length(simFol)) {
      
      # initiate new line
      newLine <- NULL
      
      # get results folder
      resFol <- paste(folder, simFol[j], sep = "/")
      print(paste("in sim folder ", j, ": ", resFol))
      
      # find stats folder and store path
      
      # get results files
      results <- list.files(resFol)
      
      # find stats folder and store path
      statsFol <- grep(pattern = c("stats"), x = results, value = T)
      
      # only csv files
      results <- grep(pattern = c(".csv"), x = results, value = T)
      
      # only the results files
      results <- grep(pattern = c(keyword), x = results, value = T)
      
      # get param values and nb of rep
      strg = resFol
      
      # strg = sub(x = strg, pattern = "*.csv", replacement = "")   # cut ".csv" out
      strg = unlist(strsplit(strg, split = "-")) # split according to "-"
      strg = strg[-c(1, 2, 3)] # take out non param elements
      
      # get param values 
      newLine <- c(newLine,
                  sub(x = strg[1], pattern = paste("predSpcf", "*", sep = ""), replacement = ""), # add to newLine everything after XparamName
                  sub(x = strg[2], pattern = paste("predOpnt", "*", sep = ""), replacement = ""),
                  sub(x = strg[3], pattern = paste("pryConvRateRatio", "*", sep = ""), replacement = ""),
                  sub(x = strg[4], pattern = paste("pryCtchProbRatio", "*", sep = ""), replacement = ""),
                  sub(x = strg[5], pattern = paste("pryOfspRatio", "*", sep = ""), replacement = ""),
                  sub(x = strg[6], pattern = paste("pryMaxConsRatio", "*", sep = ""), replacement = ""))
      
      # get nb of replicates
      newLine <- c(newLine, length(results))
      
      # initiate extinctions counters
      pry1ext = 0
      pry2ext = 0
      predExt = 0
      
      # loop over the results files
      for (m in 1:length(results)) {
        
        # print(paste("reading file ", results[m]))
        # read results
        res <- read.csv(paste(resFol, results[m], sep = "/"))
        
        # count extinctions
        if (res[dim(res)[1], 4] == 0) {pry1ext = pry1ext+1}
        if (res[dim(res)[1], 5] == 0) {pry2ext = pry2ext+1}
        if (res[dim(res)[1], 8] == 0) {predExt = predExt+1}
        
      } # end loop over files (replicates)
      
      # update new line
      newLine <- c(newLine, pry1ext/length(results), pry2ext/length(results), predExt/length(results))
      
      # browse stats folder and read stats
      stats <- list.files(paste(resFol, statsFol, sep = "/"))
      stats <- grep(pattern = c("stats"), x = stats, value = T)
      stats <- grep(pattern = c(".csv"), x = stats, value = T)
      
      res <- read.csv(paste(resFol, statsFol, stats, sep = "/"))
      
      # subset before pred introduction
      sub <- subset(res, subset = res$timeStep >= before[1] & res$timeStep < before[2])
      # get before measures
      newLine <- c(newLine, 
                   mean(sub$prey1PopulationSizeMean), max(sub$prey1PopulationSizeMean), min(sub$prey1PopulationSizeMean),
                   mean(sub$prey2PopulationSizeMean), max(sub$prey2PopulationSizeMean), min(sub$prey2PopulationSizeMean),
                   mean(sub$prey1growthRateMean/100), max(sub$prey1growthRateMean/100), min(sub$prey1growthRateMean/100),
                   mean(sub$prey2growthRateMean/100), max(sub$prey2growthRateMean/100), min(sub$prey2growthRateMean/100))
      
      # subset after
      sub <- subset(res, subset = res$timeStep >= after[1] & res$timeStep <= after[2])
      # get before measures
      newLine <- c(newLine, 
                   mean(sub$prey1PopulationSizeMean), max(sub$prey1PopulationSizeMean), min(sub$prey1PopulationSizeMean),
                   mean(sub$prey2PopulationSizeMean), max(sub$prey2PopulationSizeMean), min(sub$prey2PopulationSizeMean),
                   mean(sub$predator1PopulationSizeMean), max(sub$predator1PopulationSizeMean), min(sub$predator1PopulationSizeMean),
                   mean(sub$prey1growthRateMean/100), max(sub$prey1growthRateMean/100), min(sub$prey1growthRateMean/100),
                   mean(sub$prey2growthRateMean/100), max(sub$prey2growthRateMean/100), min(sub$prey2growthRateMean/100),
                   mean(sub$predatorGrowthRateMean/100), max(sub$predatorGrowthRateMean/100), min(sub$predatorGrowthRateMean/100),
                   mean(sub$prey1catchesMean), max(sub$prey1catchesMean), min(sub$prey1catchesMean),
                   mean(sub$prey2catchesMean), max(sub$prey2catchesMean), min(sub$prey2catchesMean))

      # rbind the line to tab
      tab <- rbind(tab, as.numeric(newLine))
      
    } # end loop over sim folders
    
    # name columns
    colnames(tab) <- headers
    
    # create local SA results folder in allStatsAndPlots and save table
    newFolderDir <- paste(folder, "allStatsAndPlots", "localSAfiles", sep = "/")
    dir.create(path = newFolderDir)
    write.csv(tab, file = paste(newFolderDir, "/stats-", content[i], ".csv", sep = ""), row.names = FALSE)
    
  } # end loop over folders
  
} # end of function

Path = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/"
Pattern = "localSA-pred"
Keyword = "Results"

mergeResults(path = Path, keyword = Keyword, pattern = Pattern)
statsResults(path = Path, keyword = Keyword, pattern = Pattern)
localSAresults(path = Path, keyword = Keyword, pattern = Pattern)

######### add missing parameter set #########

missingSet <- read.csv(file = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-PredSpecific/allStatsAndPlots/localSAfiles/stats-folder-localSA-PredSpecific.csv")
missingLine <- missingSet[2,]

# catchProba
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-CatchProb/allStatsAndPlots/localSAfiles/stats-folder-localSA-CatchProb.csv"
tab <- read.csv(file = filePath)
tab <- rbind(tab, missingLine)
write.csv(tab, file = filePath, row.names = FALSE)

# convRate
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-ConvRate/allStatsAndPlots/localSAfiles/stats-folder-localSA-ConvRate.csv"
tab <- read.csv(file = filePath)
tab <- rbind(tab, missingLine)
write.csv(tab, file = filePath, row.names = FALSE)

# maxCons
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-maxCons/allStatsAndPlots/localSAfiles/stats-folder-localSA-maxCons.csv"
tab <- read.csv(file = filePath)
tab <- rbind(tab, missingLine)
write.csv(tab, file = filePath, row.names = FALSE)

# offspring Nb
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-offspringAvg/allStatsAndPlots/localSAfiles/stats-folder-localSA-offspringAvg.csv"
tab <- read.csv(file = filePath)
tab <- rbind(tab, missingLine)
write.csv(tab, file = filePath, row.names = FALSE)

######### figures ######### 

######### Regime ########

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-PredSpecific/allStatsAndPlots/localSAfiles/"
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-PredSpecific/allStatsAndPlots/localSAfiles/stats-folder-localSA-PredSpecific.csv"
data <- read.csv(filePath)

#### Extinction Proba ####

x = as.factor(data$predSpecific)
y1 = data$prey1extFreq
y2 = data$prey2extFreq
y3 = data$prey1extFreq.1
# y1min = data$prey1PopulationSizeICinf
# y2min = data$prey2PopulationSizeICinf
# y3min = data$predator1PopulationSizeICinf
# y1max = data$prey1PopulationSizeICsup
# y2max = data$prey2PopulationSizeICsup
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  ylim(0, 1) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  # geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
  # geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_line(aes(y = y1), color = y1c) +
  # geom_line(aes(y = y2), color = y2c) +
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_point(aes(y = y1), size = 2.5, shape = 21, fill = "white", color = y1c, position = position_nudge(x = -0.15)) +
  geom_point(aes(y = y2), size = 2.5, shape = 22, fill = "white", color = y2c) +
  geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Predator Specific", y = "Extinction frequency") # +
  # scale_colour_manual(name='Populations',
  #                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
  #                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-Regime-extFreq.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### prey 2 to prey 1 density deviation (before/after) ####

y1 = data$prey2densBeforeMean/data$prey1densBeforeMean - 1
y2 = data$prey2densAfterMean/data$prey1densAfterMean - 1
# y3 = data$prey1extFreq.1
y1min = data$prey2densBeforeMin/data$prey1densBeforeMin - 1
y2min = data$prey2densAfterMin/data$prey1densAfterMin - 1
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey2densBeforeMax/data$prey1densBeforeMax - 1
y2max = data$prey2densAfterMax/data$prey1densAfterMax - 1
# y3max = data$predator1PopulationSizeICsup
y1c = "cyan"
y2c = "violet"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x, y=y1)) + 
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 1, col = y1c, position = position_nudge(x = -0.1)) +
  geom_pointrange(aes(ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 1, col = y2c, position = position_nudge(x = 0.1)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_line(aes(y = y1), color = y1c) +
  # geom_line(aes(y = y2), color = y2c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Predator Specific", y = "Average prey 2 to prey 1 density deviation") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-Regime-prey2toPrey1densDev.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### catch rate ####

y1 = data$prey1catchesMean
y2 = data$prey2catchesMean
# y3 = data$prey1extFreq.1
y1min = data$prey1catchesMin
y2min = data$prey2catchesMin
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey1catchesMax
y2max = data$prey2catchesMax
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x, y=y1)) + 
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  # geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.5, col = y1c, position = position_nudge(x = -0.1)) +
  geom_pointrange(aes(ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.5, col = y2c, position = position_nudge(x = 0.1)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_line(aes(y = y1), color = y1c) +
  # geom_line(aes(y = y2), color = y2c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Predator Specific", y = "Average catches per time step") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-Regime-catchRate.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

######### catchProba ######### 

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-CatchProb/allStatsAndPlots/localSAfiles/"
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-CatchProb/allStatsAndPlots/localSAfiles/stats-folder-localSA-CatchProb.csv"
data <- read.csv(filePath)
data <- subset(data, data$predSpecific == 0)

#### Extinction Proba ####

x = data$catchProbaRatio
y1 = data$prey1extFreq
y2 = data$prey2extFreq
y3 = data$prey1extFreq.1
# y1min = data$prey1PopulationSizeICinf
# y2min = data$prey2PopulationSizeICinf
# y3min = data$predator1PopulationSizeICinf
# y1max = data$prey1PopulationSizeICsup
# y2max = data$prey2PopulationSizeICsup
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  ylim(0, 1) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  # geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
  # geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.1)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0)) +
  geom_line(aes(y = y3), color = y3c, alpha = 0.2, position = position_nudge(x = 0.1)) +
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_point(aes(y = y1), size = 2, shape = 21, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  geom_point(aes(y = y2), size = 2, shape = 22, fill = "white", color = y2c) +
  geom_point(aes(y = y3), size = 2, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.1)) +
  labs(x = "Prey2 to prey1 catch probability ratio", y = "Extinction frequency") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-catchProba-extFreq.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### prey 2 to prey 1 density deviation (before/after) ####

y1 = data$prey2densBeforeMean/data$prey1densBeforeMean - 1
y2 = data$prey2densAfterMean/data$prey1densAfterMean - 1
# y3 = data$prey1extFreq.1
y1min = data$prey2densBeforeMin/data$prey1densBeforeMin - 1
y2min = data$prey2densAfterMin/data$prey1densAfterMin - 1
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey2densBeforeMax/data$prey1densBeforeMax - 1
y2max = data$prey2densAfterMax/data$prey1densAfterMax - 1
# y3max = data$predator1PopulationSizeICsup
y1c = "cyan"
y2c = "violet"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.1)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.1)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.8, col = y1c, position = position_nudge(x = -0.1)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.8, col = y2c, position = position_nudge(x = 0.1)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 catch probability ratio", y = "Average prey 2 to prey 1 density deviation") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-catchProba-prey2toPrey1densDev.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### catch rate ####

y1 = data$prey1catchesMean
y2 = data$prey2catchesMean
# y3 = data$prey1extFreq.1
y1min = data$prey1catchesMin
y2min = data$prey2catchesMin
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey1catchesMax
y2max = data$prey2catchesMax
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  geom_hline(yintercept = y1[3], alpha = 0.5, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y2[3], alpha = 0.5, color = "darkblue", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.05)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.1)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 21, fill = "white", size = 0.5, col = y1c, position = position_nudge(x = -0.1)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 22, fill = "white", size = 0.5, col = y2c, position = position_nudge(x = 0.1)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 catch probability ratio", y = "Average catches per time step") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-catchProba-catchRate.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

######### resources per Catch ratio ######### 

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-ConvRate/allStatsAndPlots/localSAfiles/"
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-ConvRate/allStatsAndPlots/localSAfiles/stats-folder-localSA-ConvRate.csv"
data <- read.csv(filePath)
data <- subset(data, data$predSpecific == 0 & data$predOportunistic == 0)

#### Extinction Proba ####

x = data$convRateRatio
y1 = data$prey1extFreq
y2 = data$prey2extFreq
y3 = data$prey1extFreq.1
# y1min = data$prey1PopulationSizeICinf
# y2min = data$prey2PopulationSizeICinf
# y3min = data$predator1PopulationSizeICinf
# y1max = data$prey1PopulationSizeICsup
# y2max = data$prey2PopulationSizeICsup
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  ylim(0, 1) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  # geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
  # geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.02)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0)) +
  geom_line(aes(y = y3), color = y3c, alpha = 0.2, position = position_nudge(x = 0.02)) +
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_point(aes(y = y1), size = 3, shape = 21, fill = "white", color = y1c, position = position_nudge(x = -0.02)) +
  geom_point(aes(y = y2), size = 3, shape = 22, fill = "white", color = y2c) +
  geom_point(aes(y = y3), size = 3, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.02)) +
  labs(x = "Prey2 to prey1 resources/catch ratio", y = "Extinction frequency") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-convRate-extFreq.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### prey 2 to prey 1 density deviation (before/after) ####

y1 = data$prey2densBeforeMean/data$prey1densBeforeMean - 1
y2 = data$prey2densAfterMean/data$prey1densAfterMean - 1
# y3 = data$prey1extFreq.1
y1min = data$prey2densBeforeMin/data$prey1densBeforeMin - 1
y2min = data$prey2densAfterMin/data$prey1densAfterMin - 1
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey2densBeforeMax/data$prey1densBeforeMax - 1
y2max = data$prey2densAfterMax/data$prey1densAfterMax - 1
# y3max = data$predator1PopulationSizeICsup
y1c = "cyan"
y2c = "violet"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.01)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.01)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.8, col = y1c, position = position_nudge(x = -0.01)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.8, col = y2c, position = position_nudge(x = 0.01)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 resources/catch ratio", y = "Average prey 2 to prey 1 density deviation") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-convRate-prey2toPrey1densDev.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### catch rate ####

y1 = data$prey1catchesMean
y2 = data$prey2catchesMean
# y3 = data$prey1extFreq.1
y1min = data$prey1catchesMin
y2min = data$prey2catchesMin
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey1catchesMax
y2max = data$prey2catchesMax
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  # geom_rect(aes(xmin = 0, xmax = 1, ymin = y1min[length(y1min)], ymax = y2min[length(y2min)]), alpha=0.5, fill = "lightgrey") +
  geom_hline(yintercept = y1[length(y1)], alpha = 0.5, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y2[length(y2)], alpha = 0.5, color = "darkblue", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.01)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.01)) +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 21, fill = "white", size = 0.5, col = y1c, position = position_nudge(x = -0.01)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 22, fill = "white", size = 0.5, col = y2c, position = position_nudge(x = 0.01)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 resource/catch ratio", y = "Average catches per time step") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-convRate-catchRate.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

######### average offspring nb ######### 

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-offspringAvg/allStatsAndPlots/localSAfiles/"
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-offspringAvg/allStatsAndPlots/localSAfiles/stats-folder-localSA-offspringAvg.csv"
data <- read.csv(filePath)
data <- subset(data, data$predSpecific == 0 & data$predOportunistic == 0)

#### Extinction Proba ####

x = data$preyOffspRatio
y1 = data$prey1extFreq
y2 = data$prey2extFreq
y3 = data$prey1extFreq.1
# y1min = data$prey1PopulationSizeICinf
# y2min = data$prey2PopulationSizeICinf
# y3min = data$predator1PopulationSizeICinf
# y1max = data$prey1PopulationSizeICsup
# y2max = data$prey2PopulationSizeICsup
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  ylim(0, 1) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  # geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
  # geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.2)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0)) +
  geom_line(aes(y = y3), color = y3c, alpha = 0.2, position = position_nudge(x = 0.2)) +
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_point(aes(y = y1), size = 3, shape = 21, fill = "white", color = y1c, position = position_nudge(x = -0.2)) +
  geom_point(aes(y = y2), size = 3, shape = 22, fill = "white", color = y2c) +
  geom_point(aes(y = y3), size = 3, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.2)) +
  labs(x = "Prey average offspring nb ratio", y = "Extinction frequency") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-offpringAvg-extFreq.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### prey 2 to prey 1 density deviation (before/after) ####

y1 = data$prey2densBeforeMean/data$prey1densBeforeMean - 1
y2 = data$prey2densAfterMean/data$prey1densAfterMean - 1
# y3 = data$prey1extFreq.1
y1min = data$prey2densBeforeMin/data$prey1densBeforeMin - 1
y2min = data$prey2densAfterMin/data$prey1densAfterMin - 1
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey2densBeforeMax/data$prey1densBeforeMax - 1
y2max = data$prey2densAfterMax/data$prey1densAfterMax - 1
# y3max = data$predator1PopulationSizeICsup
y1c = "cyan"
y2c = "violet"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.1)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.1)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.8, col = y1c, position = position_nudge(x = -0.1)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.8, col = y2c, position = position_nudge(x = 0.1)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 average offspring ratio", y = "Average prey 2 to prey 1 density deviation") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-offsprinAvg-prey2toPrey1densDev.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### catch rate ####

y1 = data$prey1catchesMean
y2 = data$prey2catchesMean
# y3 = data$prey1extFreq.1
y1min = data$prey1catchesMin
y2min = data$prey2catchesMin
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey1catchesMax
y2max = data$prey2catchesMax
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  # geom_rect(aes(xmin = 0, xmax = 1, ymin = y1min[length(y1min)], ymax = y2min[length(y2min)]), alpha=0.5, fill = "lightgrey") +
  geom_hline(yintercept = y1[1], alpha = 0.5, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y2[1], alpha = 0.5, color = "darkblue", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.1)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.1)) +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 21, fill = "white", size = 0.5, col = y1c, position = position_nudge(x = -0.1)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 22, fill = "white", size = 0.5, col = y2c, position = position_nudge(x = 0.1)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 average offspring nb ratio", y = "Average catches per time step") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-offspringAvg-catchRate.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

######### maximum consumption ######### 

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-maxCons/allStatsAndPlots/localSAfiles/"
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-maxCons/allStatsAndPlots/localSAfiles/stats-folder-localSA-maxCons.csv"
data <- read.csv(filePath)
data <- subset(data, data$predSpecific == 0 & data$predOportunistic == 0)

#### Extinction Proba ####

x = data$maxConsRatio
y1 = data$prey1extFreq
y2 = data$prey2extFreq
y3 = data$prey1extFreq.1
# y1min = data$prey1PopulationSizeICinf
# y2min = data$prey2PopulationSizeICinf
# y3min = data$predator1PopulationSizeICinf
# y1max = data$prey1PopulationSizeICsup
# y2max = data$prey2PopulationSizeICsup
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  ylim(0, 1) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  # geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
  # geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.05)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0)) +
  geom_line(aes(y = y3), color = y3c, alpha = 0.2, position = position_nudge(x = 0.05)) +
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_point(aes(y = y1), size = 3, shape = 21, fill = "white", color = y1c, position = position_nudge(x = -0.05)) +
  geom_point(aes(y = y2), size = 3, shape = 22, fill = "white", color = y2c) +
  geom_point(aes(y = y3), size = 3, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.05)) +
  labs(x = "Prey2 to prey1 max consumption ratio", y = "Extinction frequency") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-maxCons-extFreq.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### prey 2 to prey 1 density deviation (before/after) ####

y1 = data$prey2densBeforeMean/data$prey1densBeforeMean - 1
y2 = data$prey2densAfterMean/data$prey1densAfterMean - 1
# y3 = data$prey1extFreq.1
y1min = data$prey2densBeforeMin/data$prey1densBeforeMin - 1
y2min = data$prey2densAfterMin/data$prey1densAfterMin - 1
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey2densBeforeMax/data$prey1densBeforeMax - 1
y2max = data$prey2densAfterMax/data$prey1densAfterMax - 1
# y3max = data$predator1PopulationSizeICsup
y1c = "cyan"
y2c = "violet"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.01)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.01)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.8, col = y1c, position = position_nudge(x = -0.01)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.8, col = y2c, position = position_nudge(x = 0.01)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 max consumption ratio", y = "Average prey 2 to prey 1 density deviation") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-maxCons-prey2toPrey1densDev.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

fig <- ggplot(data[1:3,], aes(x[1:3])) + 
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  # geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_line(aes(y = y1[1:3]), color = y1c, alpha = 0.2, position = position_nudge(x = -0.02)) +
  geom_line(aes(y = y2[1:3]), color = y2c, alpha = 0.2, position = position_nudge(x = 0.02)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1[1:3], ymin = y1min[1:3], ymax = y1max[1:3]), shape = 24, fill = "white", size = 0.8, col = y1c, position = position_nudge(x = -0.02)) +
  geom_pointrange(aes(y = y2[1:3], ymin = y2min[1:3], ymax = y2max[1:3]), shape = 25, fill = "white", size = 0.8, col = y2c, position = position_nudge(x = 0.02)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 max consumption ratio", y = "Average prey 2 to prey 1 density deviation") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-maxCons-prey2toPrey1densDev-under1.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

fig <- ggplot(data[3:6,], aes(x[3:6])) + 
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_line(aes(y = y1[3:6]), color = y1c, alpha = 0.2, position = position_nudge(x = -0.02)) +
  geom_line(aes(y = y2[3:6]), color = y2c, alpha = 0.2, position = position_nudge(x = 0.02)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1[3:6], ymin = y1min[3:6], ymax = y1max[3:6]), shape = 24, fill = "white", size = 0.8, col = y1c, position = position_nudge(x = -0.02)) +
  geom_pointrange(aes(y = y2[3:6], ymin = y2min[3:6], ymax = y2max[3:6]), shape = 25, fill = "white", size = 0.8, col = y2c, position = position_nudge(x = 0.02)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 max consumption ratio", y = "Average prey 2 to prey 1 density deviation") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-maxCons-prey2toPrey1densDev-over1.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### catch rate ####

y1 = data$prey1catchesMean
y2 = data$prey2catchesMean
# y3 = data$prey1extFreq.1
y1min = data$prey1catchesMin
y2min = data$prey2catchesMin
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey1catchesMax
y2max = data$prey2catchesMax
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
# y3c = "orange"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  # geom_rect(aes(xmin = 0, xmax = 1, ymin = y1min[length(y1min)], ymax = y2min[length(y2min)]), alpha=0.5, fill = "lightgrey") +
  geom_hline(yintercept = y1[3], alpha = 0.5, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y2[3], alpha = 0.5, color = "darkblue", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.01)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.01)) +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 21, fill = "white", size = 0.5, col = y1c, position = position_nudge(x = -0.01)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 22, fill = "white", size = 0.5, col = y2c, position = position_nudge(x = 0.01)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 max consumption ratio", y = "Average catches per time step") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-maxCons-catchRate.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

#### New figures for density ####

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-offspringAvg/allStatsAndPlots/localSAfiles/"
filePath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-offspringAvg/allStatsAndPlots/localSAfiles/stats-folder-localSA-offspringAvg.csv"
data <- read.csv(filePath)
data <- subset(data, data$predSpecific == 0 & data$predOportunistic == 0)

x = data$preyOffspRatio
y1 = data$prey1densBeforeMean
y2 = data$prey1densAfterMean
y3 = data$prey2densBeforeMean
y4 = data$prey2densAfterMean
y1min = data$prey1densBeforeMin
y2min = data$prey1densAfterMin
y3min = data$prey2densBeforeMin
y4min = data$prey2densAfterMin
y1max = data$prey1densBeforeMax
y2max = data$prey1densAfterMax
y3max = data$prey2densBeforeMax
y4max = data$prey2densAfterMax
y1c = "pink"
y2c = "red"
y3c = "cyan"
y4c = "blue"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  # geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y1[1], alpha = 0.5, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y3[1], alpha = 0.5, color = "darkblue", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.15)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.05)) +
  geom_line(aes(y = y3), color = y3c, alpha = 0.2, position = position_nudge(x = -0.05)) +
  geom_line(aes(y = y4), color = y4c, alpha = 0.2, position = position_nudge(x = 0.15)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.8, col = y1c, position = position_nudge(x = -0.15)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.8, col = y2c, position = position_nudge(x = 0.05)) +
  geom_pointrange(aes(y = y3, ymin = y3min, ymax = y3max), shape = 24, fill = "white", size = 0.8, col = y3c, position = position_nudge(x = -0.05)) +
  geom_pointrange(aes(y = y4, ymin = y4min, ymax = y4max), shape = 25, fill = "white", size = 0.8, col = y4c, position = position_nudge(x = 0.15)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 average offspring ratio", y = "Average prey 2 to prey 1 density") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-offsprinAvg-prey2toPrey1densDev-new.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-CatchProb/allStatsAndPlots/localSAfiles/"
filePath = paste(folderPath, "stats-folder-localSA-CatchProb.csv", sep = "")
data <- read.csv(filePath)
data <- subset(data, data$predSpecific == 0 & data$predOportunistic == 0)

x = data$catchProbaRatio
y1 = data$prey1densBeforeMean
y2 = data$prey1densAfterMean
y3 = data$prey2densBeforeMean
y4 = data$prey2densAfterMean
y1min = data$prey1densBeforeMin
y2min = data$prey1densAfterMin
y3min = data$prey2densBeforeMin
y4min = data$prey2densAfterMin
y1max = data$prey1densBeforeMax
y2max = data$prey1densAfterMax
y3max = data$prey2densBeforeMax
y4max = data$prey2densAfterMax
y1c = "pink"
y2c = "red"
y3c = "cyan"
y4c = "blue"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  # geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y1[1], alpha = 0.5, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y3[1], alpha = 0.5, color = "darkblue", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.15)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.05)) +
  geom_line(aes(y = y3), color = y3c, alpha = 0.2, position = position_nudge(x = -0.05)) +
  geom_line(aes(y = y4), color = y4c, alpha = 0.2, position = position_nudge(x = 0.15)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.8, col = y1c, position = position_nudge(x = -0.15)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.8, col = y2c, position = position_nudge(x = 0.05)) +
  geom_pointrange(aes(y = y3, ymin = y3min, ymax = y3max), shape = 24, fill = "white", size = 0.8, col = y3c, position = position_nudge(x = -0.05)) +
  geom_pointrange(aes(y = y4, ymin = y4min, ymax = y4max), shape = 25, fill = "white", size = 0.8, col = y4c, position = position_nudge(x = 0.15)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 catch probability ratio", y = "Average prey 2 to prey 1 density") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-catchProba-prey2toPrey1densDev-new.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-ConvRate/allStatsAndPlots/localSAfiles/"
filePath = paste(folderPath, "stats-folder-localSA-ConvRate.csv", sep = "")
data <- read.csv(filePath)
data <- subset(data, data$predSpecific == 0 & data$predOportunistic == 0)

x = data$convRateRatio
y1 = data$prey1densBeforeMean
y2 = data$prey1densAfterMean
y3 = data$prey2densBeforeMean
y4 = data$prey2densAfterMean
y1min = data$prey1densBeforeMin
y2min = data$prey1densAfterMin
y3min = data$prey2densBeforeMin
y4min = data$prey2densAfterMin
y1max = data$prey1densBeforeMax
y2max = data$prey1densAfterMax
y3max = data$prey2densBeforeMax
y4max = data$prey2densAfterMax
y1c = "pink"
y2c = "red"
y3c = "cyan"
y4c = "blue"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  # geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y1[1], alpha = 0.5, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y3[1], alpha = 0.5, color = "darkblue", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.03)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.01)) +
  geom_line(aes(y = y3), color = y3c, alpha = 0.2, position = position_nudge(x = -0.01)) +
  geom_line(aes(y = y4), color = y4c, alpha = 0.2, position = position_nudge(x = 0.03)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.5, col = y1c, position = position_nudge(x = -0.03)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.5, col = y2c, position = position_nudge(x = 0.01)) +
  geom_pointrange(aes(y = y3, ymin = y3min, ymax = y3max), shape = 24, fill = "white", size = 0.5, col = y3c, position = position_nudge(x = -0.01)) +
  geom_pointrange(aes(y = y4, ymin = y4min, ymax = y4max), shape = 25, fill = "white", size = 0.5, col = y4c, position = position_nudge(x = 0.03)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 average resource/catch ratio", y = "Average prey 2 to prey 1 density") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-convRate-prey2toPrey1densDev-new.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-maxCons/allStatsAndPlots/localSAfiles/"
filePath = paste(folderPath, "stats-folder-localSA-maxCons.csv", sep = "")
data <- read.csv(filePath)
data <- subset(data, data$predSpecific == 0 & data$predOportunistic == 0)

x = data$maxConsRatio
y1 = data$prey1densBeforeMean
y2 = data$prey1densAfterMean
y3 = data$prey2densBeforeMean
y4 = data$prey2densAfterMean
y1min = data$prey1densBeforeMin
y2min = data$prey1densAfterMin
y3min = data$prey2densBeforeMin
y4min = data$prey2densAfterMin
y1max = data$prey1densBeforeMax
y2max = data$prey1densAfterMax
y3max = data$prey2densBeforeMax
y4max = data$prey2densAfterMax
y1c = "pink"
y2c = "red"
y3c = "cyan"
y4c = "blue"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  # geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y1[3], alpha = 0.5, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y3[3], alpha = 0.5, color = "darkblue", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.04)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.01)) +
  geom_line(aes(y = y3), color = y3c, alpha = 0.2, position = position_nudge(x = -0.01)) +
  geom_line(aes(y = y4), color = y4c, alpha = 0.2, position = position_nudge(x = 0.04)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.5, col = y1c, position = position_nudge(x = -0.04)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.5, col = y2c, position = position_nudge(x = 0.01)) +
  geom_pointrange(aes(y = y3, ymin = y3min, ymax = y3max), shape = 24, fill = "white", size = 0.5, col = y3c, position = position_nudge(x = -0.01)) +
  geom_pointrange(aes(y = y4, ymin = y4min, ymax = y4max), shape = 25, fill = "white", size = 0.5, col = y4c, position = position_nudge(x = 0.04)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Prey2 to prey1 max consumption ratio", y = "Average prey 2 to prey 1 density") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-maxCons-prey2toPrey1densDev-new.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)

folderPath = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/localSA/folder-localSA-PredSpecific/allStatsAndPlots/localSAfiles/"
filePath = paste(folderPath, "stats-folder-localSA-PredSpecific.csv", sep = "")
data <- read.csv(filePath)
# data <- subset(data, data$predSpecific == 0 & data$predOportunistic == 0)

x = as.factor(data$predSpecific)
y1 = data$prey1densBeforeMean
y2 = data$prey1densAfterMean
y3 = data$prey2densBeforeMean
y4 = data$prey2densAfterMean
y1min = data$prey1densBeforeMin
y2min = data$prey1densAfterMin
y3min = data$prey2densBeforeMin
y4min = data$prey2densAfterMin
y1max = data$prey1densBeforeMax
y2max = data$prey1densAfterMax
y3max = data$prey2densBeforeMax
y4max = data$prey2densAfterMax
y1c = "pink"
y2c = "red"
y3c = "cyan"
y4c = "blue"
# tIntro = 210

fig <- ggplot(data, aes(x)) + 
  # geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  # geom_hline(yintercept = -1, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y1[3], alpha = 0.5, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = y3[3], alpha = 0.5, color = "darkblue", linetype = "dashed") +
  geom_line(aes(y = y1), color = y1c, alpha = 0.2, position = position_nudge(x = -0.075)) +
  geom_line(aes(y = y2), color = y2c, alpha = 0.2, position = position_nudge(x = 0.025)) +
  geom_line(aes(y = y3), color = y3c, alpha = 0.2, position = position_nudge(x = -0.025)) +
  geom_line(aes(y = y4), color = y4c, alpha = 0.2, position = position_nudge(x = 0.075)) +
  # geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(data$prey2PopulationSizeMean)), alpha=0.5, fill = "lightgrey") +
  geom_pointrange(aes(y = y1, ymin = y1min, ymax = y1max), shape = 24, fill = "white", size = 0.5, col = y1c, position = position_nudge(x = -0.075)) +
  geom_pointrange(aes(y = y2, ymin = y2min, ymax = y2max), shape = 25, fill = "white", size = 0.5, col = y2c, position = position_nudge(x = 0.025)) +
  geom_pointrange(aes(y = y3, ymin = y3min, ymax = y3max), shape = 24, fill = "white", size = 0.5, col = y3c, position = position_nudge(x = -0.025)) +
  geom_pointrange(aes(y = y4, ymin = y4min, ymax = y4max), shape = 25, fill = "white", size = 0.5, col = y4c, position = position_nudge(x = 0.075)) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 24, fill = "white", color = y1c, position = position_nudge(x = -0.1)) +
  # geom_point(aes(y = y2), size = 2.5, shape = 25, fill = "white", color = y2c, position = position_nudge(x = 0.1)) # +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c, position = position_nudge(x = 0.15)) +
  labs(x = "Predators specific", y = "Average prey 2 to prey 1 density") # +
# scale_colour_manual(name='Populations',
#                     breaks=c('Prey 1', 'Prey 2', 'Predator'),
#                     values=c(y1c, y2c, y3c))

# save plot in this folder
ggsave(filename = "localSA-predRegime-prey2toPrey1densDev-new.pdf", path = folderPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
