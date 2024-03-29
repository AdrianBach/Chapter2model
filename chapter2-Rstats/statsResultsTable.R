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
    folder = paste(path, content[i], sep = "")
    
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
  
  # create a global stats folder
  dir.create(path = paste(path, "allStatsAndPlots", sep = "/"))
  
  # create a Results folder
  dir.create(path = paste(path, "/allStatsAndPlots/", keyword, "Files", sep = ""))
  
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

    # and copy in the global stats folder
    file.copy(from = paste(folder, "/stats", keyword, "-", content[i], ".csv", sep = ""), to = paste(path, "/allStatsAndPlots/", keyword, "Files", sep = ""))
    
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
    ggsave(filename = paste("stats", keyword, "-", content[i], ".pdf", sep = ""), path = folder, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
    # and copy in the global stats folder
    file.copy(from = paste(folder, "/stats", keyword, "-", content[i], ".pdf", sep = ""), to = paste(path, "/allStatsAndPlots/", keyword, "Files", sep = ""))
    
    } # end of loop over folders
  } # end of function

plotDensity <- function(path, Pattern) {
  
  # create a folder for figures
  dirPath = paste(path, "density", sep = "")
  dir.create(path = dirPath)
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c(Pattern), x = content, value = T)
  
  # only the csv
  content <- grep(pattern = c(".csv"), x = content, value = T)
  
  # loop over content
  for (i in 1:length(content)) {
    
    # # path to folder
    # # folder = paste("~/", path, content[i], "/stats-", content[i], sep = "")
    # folder = paste(path, "/", content[i], "/stats-", content[i], sep = "")
    
    # # get content
    # results = list.files(folder)
    
    # # select only results csv files
    # results <- grep(pattern = c(keyword), x = results, value = T)
    # results <- grep(pattern = c("stats"), x = results, value = T)
    # results <- grep(pattern = c(".csv"), x = results, value = T)
    
    # # copy in the global stats folder
    # file.copy(from = paste(folder, results, sep = "/"), to = paste(path, "allStatsAndPlots", sep = "/"))
  
    # plot and save figure
    data <- read.csv(paste(path, content[i], sep = ""))
    
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
      geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(y1max)), alpha=0.5, fill = "lightgrey") +
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
    
    # # save plot in this folder
    # ggsave(filename = paste("stats", keyword, "-", content[i], ".pdf", sep = ""), path = folder, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
    # and in the global stats folder
    ggsave(filename = paste("density", "-", content[i], ".pdf", sep = ""), path = dirPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
  } # end of loop over content
} # end of function

plotCatches <- function(path, Pattern) {
  
  # create a folder for figures
  dirPath = paste(path, "catches", sep = "")
  dir.create(path = dirPath)
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c(Pattern), x = content, value = T)
  
  # only the csv
  content <- grep(pattern = c(".csv"), x = content, value = T)
  
  # loop over content
  for (i in 1:length(content)) {
    
    # # path to folder
    # # folder = paste("~/", path, content[i], "/stats-", content[i], sep = "")
    # folder = paste(path, "/", content[i], "/stats-", content[i], sep = "")
    
    # # get content
    # results = list.files(folder)
    
    # # select only results csv files
    # results <- grep(pattern = c(keyword), x = results, value = T)
    # results <- grep(pattern = c("stats"), x = results, value = T)
    # results <- grep(pattern = c(".csv"), x = results, value = T)
    
    # # copy in the global stats folder
    # file.copy(from = paste(folder, results, sep = "/"), to = paste(path, "allStatsAndPlots", sep = "/"))
    
    # plot and save figure
    data <- read.csv(paste(path, content[i], sep = ""))
    
    x = data$timeStep
    y1 = data$prey1catchesMean
    y2 = data$prey2catchesMean
    y1min = data$prey1catchesICinf
    y2min = data$prey2catchesICinf
    y1max = data$prey1catchesICsup
    y2max = data$prey2catchesICsup
    y1c = "red"
    y2c = "blue"
    tIntro = 210
    
    fig <- ggplot(data, aes(x)) + 
      geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(y1max)), alpha=0.5, fill = "lightgrey") +
      geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
      geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
      # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
      geom_line(aes(y = y1), color = y1c) +
      geom_line(aes(y = y2), color = y2c) +
      # geom_line(aes(y = y3), color = y3c) +
      # geom_point(aes(y = y1), size = 2.5, shape = 21, fill = "white", color = y1c) +
      # geom_point(aes(y = y2), size = 2.5, shape = 22, fill = "white", color = y2c) +
      # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c) +
      labs(x = "Time steps", y = "Nb of catches") # +
      # scale_colour_manual(name='Populations',
                          # breaks=c('Prey 1', 'Prey 2', 'Predator'),
                          # values=c(y1c, y2c, y3c))
    
    # # save plot in this folder
    # ggsave(filename = paste("stats", keyword, "-", content[i], ".pdf", sep = ""), path = folder, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
    # and in the global stats folder
    ggsave(filename = paste("catches", "-", content[i], ".pdf", sep = ""), path = dirPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
  } # end of loop over content
} # end of function

plotGR <- function(path, Pattern) {
  
  # create a folder for figures
  dirPath = paste(path, "growthRates", sep = "")
  dir.create(path = dirPath)
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c(Pattern), x = content, value = T)
  
  # only the csv
  content <- grep(pattern = c(".csv"), x = content, value = T)
  
  # loop over content
  for (i in 1:length(content)) {
    
    # # path to folder
    # # folder = paste("~/", path, content[i], "/stats-", content[i], sep = "")
    # folder = paste(path, "/", content[i], "/stats-", content[i], sep = "")
    
    # # get content
    # results = list.files(folder)
    
    # # select only results csv files
    # results <- grep(pattern = c(keyword), x = results, value = T)
    # results <- grep(pattern = c("stats"), x = results, value = T)
    # results <- grep(pattern = c(".csv"), x = results, value = T)
    
    # # copy in the global stats folder
    # file.copy(from = paste(folder, results, sep = "/"), to = paste(path, "allStatsAndPlots", sep = "/"))
    
    # plot and save figure
    data <- read.csv(paste(path, content[i], sep = ""))
    
    x = data$timeStep
    y1 = data$prey1growthRateMean/100
    y2 = data$prey2growthRateMean/100
    y3 = data$predatorGrowthRateMean/100
    y1min = data$prey1growthRateICinf/100
    y2min = data$prey2growthRateICinf/100
    y3min = data$predatorGrowthRateICinf/100
    y1max = data$prey1growthRateICsup/100
    y2max = data$prey2growthRateICsup/100
    y3max = data$predatorGrowthRateICsup/100
    y1c = "red"
    y2c = "blue"
    y3c = "orange"
    tIntro = 210
    
    fig <- ggplot(data, aes(x)) + 
      geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 1.05*min(y1min), ymax = 1.05*max(y1max)), alpha=0.5, fill = "lightgrey") +
      geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
      geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
      geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
      geom_line(aes(y = y3), color = y3c) +
      geom_line(aes(y = y1), color = y1c) +
      geom_line(aes(y = y2), color = y2c) +
      # geom_point(aes(y = y1), size = 2.5, shape = 21, fill = "white", color = y1c) +
      # geom_point(aes(y = y2), size = 2.5, shape = 22, fill = "white", color = y2c) +
      # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c) +
      labs(x = "Time steps", y = "Growth rate") # +
    # scale_colour_manual(name='Populations',
    # breaks=c('Prey 1', 'Prey 2', 'Predator'),
    # values=c(y1c, y2c, y3c))
    
    # # save plot in this folder
    # ggsave(filename = paste("stats", keyword, "-", content[i], ".pdf", sep = ""), path = folder, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
    # and in the global stats folder
    ggsave(filename = paste("growthRate", "-", content[i], ".pdf", sep = ""), path = dirPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
  } # end of loop over content
} # end of function

plotRelativeGR <- function(path, Pattern) {
  
  # create a folder for figures
  dirPath = paste(path, "relativeGrowthRates", sep = "")
  dir.create(path = dirPath)
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c(Pattern), x = content, value = T)
  
  # only the csv
  content <- grep(pattern = c(".csv"), x = content, value = T)
  
  # loop over content
  for (i in 1:length(content)) {
    
    # # path to folder
    # # folder = paste("~/", path, content[i], "/stats-", content[i], sep = "")
    # folder = paste(path, "/", content[i], "/stats-", content[i], sep = "")
    
    # # get content
    # results = list.files(folder)
    
    # # select only results csv files
    # results <- grep(pattern = c(keyword), x = results, value = T)
    # results <- grep(pattern = c("stats"), x = results, value = T)
    # results <- grep(pattern = c(".csv"), x = results, value = T)
    
    # # copy in the global stats folder
    # file.copy(from = paste(folder, results, sep = "/"), to = paste(path, "allStatsAndPlots", sep = "/"))
    
    # plot and save figure
    data <- read.csv(paste(path, content[i], sep = ""))
    len <- dim(data)[1]
    data <- data[-1,]
    data <- data[-len,]
    
    x = data$timeStep
    y1 = (data$prey1growthRateMean/100)/(data$prey2growthRateMean/100)
    # y2 = data$prey2growthRateMean/100
    y3 = data$predatorGrowthRateMean/100
    y1min = (data$prey1growthRateICinf/100)/(data$prey2growthRateICinf/100)
    # y2min = data$prey2growthRateICinf/100
    y3min = data$predatorGrowthRateICinf/100
    y1max = (data$prey1growthRateICsup/100)/(data$prey2growthRateICsup/100)
    # y2max = data$prey2growthRateICsup/100
    y3max = data$predatorGrowthRateICsup/100
    y1c = "black"
    # y2c = "blue"
    y3c = "orange"
    tIntro = 210
    
    fig <- ggplot(data, aes(x)) + 
      geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 1.05*min(1/y1), ymax = 1.05*max(1/y1)), alpha=0.5, fill = "lightgrey") +
      # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
      # geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
      # geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
      # geom_line(aes(y = y3), color = y3c) +
      geom_line(aes(y = 1), color = "darkgreen", linetype="dashed", size = 1) +
      geom_line(aes(y = 1/y1), color = y1c) +
      # geom_point(aes(y = y1), size = 2.5, shape = 21, fill = "white", color = y1c) +
      # geom_point(aes(y = y2), size = 2.5, shape = 22, fill = "white", color = y2c) +
      # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c) +
      labs(x = "Time steps", y = "Prey1 to prey2 relative mean growth rate") # +
    # scale_colour_manual(name='Populations',
    # breaks=c('Prey 1', 'Prey 2', 'Predator'),
    # values=c(y1c, y2c, y3c))
    
    # # save plot in this folder
    # ggsave(filename = paste("stats", keyword, "-", content[i], ".pdf", sep = ""), path = folder, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
    # and in the global stats folder
    ggsave(filename = paste("relativeGR", "-", content[i], ".pdf", sep = ""), path = dirPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
  } # end of loop over content
} # end of function

plotRelativeDensities <- function(path, Pattern) {
  
  # create a folder for figures
  dirPath = paste(path, "relativeDensities", sep = "")
  dir.create(path = dirPath)
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c(Pattern), x = content, value = T)
  
  # only the csv
  content <- grep(pattern = c(".csv"), x = content, value = T)
  
  # loop over content
  for (i in 1:length(content)) {
    
    # # path to folder
    # # folder = paste("~/", path, content[i], "/stats-", content[i], sep = "")
    # folder = paste(path, "/", content[i], "/stats-", content[i], sep = "")
    
    # # get content
    # results = list.files(folder)
    
    # # select only results csv files
    # results <- grep(pattern = c(keyword), x = results, value = T)
    # results <- grep(pattern = c("stats"), x = results, value = T)
    # results <- grep(pattern = c(".csv"), x = results, value = T)
    
    # # copy in the global stats folder
    # file.copy(from = paste(folder, results, sep = "/"), to = paste(path, "allStatsAndPlots", sep = "/"))
    
    # plot and save figure
    data <- read.csv(paste(path, content[i], sep = ""))
    
    x = data$timeStep
    y1 = (data$prey1PopulationSizeMean)/(data$prey2PopulationSizeMean)
    # y2 = data$prey2PopulationSizeMean
    y3 = data$predatorPopulationSizeMean
    y1min = (data$prey1PopulationSizeICinf)/(data$prey2PopulationSizeICinf)
    # y2min = data$prey2PopulationSizeICinf
    y3min = data$predatorPopulationSizeICinf
    y1max = (data$prey1PopulationSizeICsup)/(data$prey2PopulationSizeICsup)
    # y2max = data$prey2PopulationSizeICsup
    y3max = data$predatorPopulationSizeICsup
    y1c = "black"
    # y2c = "blue"
    y3c = "orange"
    tIntro = 210
    
    fig <- ggplot(data, aes(x)) + 
      geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0.99*min(1/y1), ymax = 1.05*max(1/y1)), alpha=0.5, fill = "lightgrey") +
      # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
      # geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
      # geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
      # geom_line(aes(y = y3), color = y3c) +
      geom_line(aes(y = 1), color = "darkgreen", linetype="dashed", size = 1) +
      geom_line(aes(y = 1/y1), color = y1c) +
      # geom_point(aes(y = y1), size = 2.5, shape = 21, fill = "white", color = y1c) +
      # geom_point(aes(y = y2), size = 2.5, shape = 22, fill = "white", color = y2c) +
      # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c) +
      labs(x = "Time steps", y = "Prey2 to prey1 relative mean densities") # +
    # scale_colour_manual(name='Populations',
    # breaks=c('Prey 1', 'Prey 2', 'Predator'),
    # values=c(y1c, y2c, y3c))
    
    # # save plot in this folder
    # ggsave(filename = paste("stats", keyword, "-", content[i], ".pdf", sep = ""), path = folder, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
    # and in the global stats folder
    ggsave(filename = paste("relativeDensity", "-", content[i], ".pdf", sep = ""), path = dirPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
  } # end of loop over content
} # end of function

plotRelativeCatches <- function(path, Pattern) {
  
  # create a folder for figures
  dirPath = paste(path, "relativeCatches", sep = "")
  dir.create(path = dirPath)
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the folders
  content <- grep(pattern = c(Pattern), x = content, value = T)
  
  # only the csv
  content <- grep(pattern = c(".csv"), x = content, value = T)
  
  # loop over content
  for (i in 1:length(content)) {
    
    # # path to folder
    # # folder = paste("~/", path, content[i], "/stats-", content[i], sep = "")
    # folder = paste(path, "/", content[i], "/stats-", content[i], sep = "")
    
    # # get content
    # results = list.files(folder)
    
    # # select only results csv files
    # results <- grep(pattern = c(keyword), x = results, value = T)
    # results <- grep(pattern = c("stats"), x = results, value = T)
    # results <- grep(pattern = c(".csv"), x = results, value = T)
    
    # # copy in the global stats folder
    # file.copy(from = paste(folder, results, sep = "/"), to = paste(path, "allStatsAndPlots", sep = "/"))
    
    # plot and save figure
    data <- read.csv(paste(path, content[i], sep = ""))
    
    x = data$timeStep
    y1 = (data$prey1catchesMean)/(data$prey2catchesMean)
    # y2 = data$prey2catchesMean
    y3 = data$predatorcatchesMean
    y1min = (data$prey1catchesICinf)/(data$prey2catchesICinf)
    # y2min = data$prey2catchesICinf
    y3min = data$predatorcatchesICinf
    y1max = (data$prey1catchesICsup)/(data$prey2catchesICsup)
    # y2max = data$prey2catchesICsup
    y3max = data$predatorcatchesICsup
    y1c = "black"
    # y2c = "blue"
    y3c = "orange"
    tIntro = 210
    
    fig <- ggplot(data, aes(x)) + 
      geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0.99*min(1/y1, na.rm = T), ymax = 1.05*max(1/y1, na.rm = T)), alpha=0.5, fill = "lightgrey") +
      # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
      # geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
      # geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
      # geom_line(aes(y = y3), color = y3c) +
      geom_line(aes(y = 1), color = "darkgreen", linetype="dashed", size = 1) +
      geom_line(aes(y = 1/y1), color = y1c) +
      # geom_point(aes(y = y1), size = 2.5, shape = 21, fill = "white", color = y1c) +
      # geom_point(aes(y = y2), size = 2.5, shape = 22, fill = "white", color = y2c) +
      # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c) +
      labs(x = "Time steps", y = "Prey2 to prey1 relative mean catches") # +
    # scale_colour_manual(name='Populations',
    # breaks=c('Prey 1', 'Prey 2', 'Predator'),
    # values=c(y1c, y2c, y3c))
    
    # # save plot in this folder
    # ggsave(filename = paste("stats", keyword, "-", content[i], ".pdf", sep = ""), path = folder, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
    # and in the global stats folder
    ggsave(filename = paste("relativeCatches", "-", content[i], ".pdf", sep = ""), path = dirPath, plot = fig, width = 6.22, height = 5.73, limitsize = TRUE)
    
  } # end of loop over content
} # end of function

gridFigDensity <- function(path, XparamName, YparamName) {
  
  # first merge all results
  
  # get the directory content
  # content <- list.files(paste("~/", path, sep = ""))
  content <- list.files(path)
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the tables
  content <- grep(pattern = ".csv", x = content, value = T)
  
  # only the ones begining with "stats"
  content <- grep(pattern = "stats", x = content, value = T)
  
  # loop over the content
  for (i in 1:length(content)) {
    
    # get the value of XparamName
    strg = content[i]
    strg = sub(x = strg, pattern = "*.csv", replacement = "")   # cut ".csv" out
    strg = unlist(strsplit(strg, split = "-")) # split according to "-"
    Xparam = strg[which(str_detect(strg, XparamName) == TRUE)] # search for XparamName in the list
    Xparam = sub(x = Xparam, pattern = paste(XparamName, "*", sep = ""), replacement = "") # get everything after XparamName
    
    # get the value of YparamName
    Yparam = strg[which(str_detect(strg, YparamName) == TRUE)] # search for XparamName in the list
    Yparam = sub(x = Yparam, pattern = paste(YparamName, "*", sep = ""), replacement = "") # get everything after XparamName
    
    # read one file to get column number and headers
    res <- read.csv(paste(path, content[i], sep = "/"))
      
    # if first round, create merge table
    if (i == 1) {
      
      # create an empty table
      headers <- colnames(res)
      headers <- c(XparamName, YparamName, colnames(res))
      tab <- data.frame(matrix(ncol = length(headers), nrow = 0))
    }
      
    # bind the replicate number to the table
    x <- rep(Xparam, times = dim(res)[1])
    y <- rep(Yparam, times = dim(res)[1])
    res <- cbind(x, y, res)
    
    tab <- rbind(tab, res) 
    
  } # end loop over folders
    
  # name columns
  colnames(tab) <- headers
  
  # save table in case
  write.csv(tab, file = paste(path,"/", XparamName, "X", YparamName, "MergedTable.csv", sep = ""), row.names = FALSE)
  
  # plot
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
    geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(y1max)), alpha=0.5, fill = "lightgrey") +
    geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
    geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
    geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
    geom_line(aes(y = y1), color = y1c) +
    geom_line(aes(y = y2), color = y2c) +
    geom_line(aes(y = y3), color = y3c)
  
  ggp <- fig + facet_grid(data$pry1cons~data$prdCtchProb1) +
    labs(x = "Time steps", y = "Population size")
  
  # save plot in this folder
  ggsave(filename = paste(XparamName, "X", YparamName, ".pdf", sep = ""), path = path, plot = ggp, width = 6.22, height = 5.73, limitsize = TRUE)
  
  
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
path = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/findEqNarrow/"
# path = getwd()
keyword = "Results"
Pattern = "findEqNarrow-size"

mergeResults(path, keyword, Pattern)
statsResults(path, keyword, Pattern)

stabPeriod = 5
stabThreshold = 10

posStabNoPred(path, keyword, Pattern, stabThreshold, stabPeriod)

path = "/home/adrian/Documents/GitKraken/Chapter2model/findEqNarrow/allStatsAndPlots/ResultsFiles/"

plotCatches(path, Pattern)
plotGR(path, Pattern)
plotRelativeGR(path, Pattern)
plotRelativeDensities(path, Pattern)
plotRelativeCatches(path, Pattern)

path = "/Users/adrianbach/Desktop/PhD/GitKraken/Chapter2model/findEqNarrow/allStatsAndPlots/ResultsFiles/"
XparamName = "pry1cons"
YparamName = "prdCtchProb1"

gridFig(path, XparamName, YparamName)
