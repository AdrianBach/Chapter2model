path = "/home/adrian/Documents/GitKraken/Chapter2model/findEqNarrow/allStatsAndPlots/ResultsFiles/"
data <- read.csv(paste(path, "pry1consXprdCtchProb1MergedTable.csv", sep = ""))

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
  geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0.95*min(1/y1), ymax = 1.05*max(1/y1)), alpha=0.5, fill = "lightgrey") +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  # geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
  # geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
  # geom_line(aes(y = y3), color = y3c) +
  geom_line(aes(y = 1), color = "darkgreen", linetype="dashed", size = 1) +
  geom_line(aes(y = 1/y1), color = y1c) +
  # geom_point(aes(y = y1), size = 2.5, shape = 21, fill = "white", color = y1c) +
  # geom_point(aes(y = y2), size = 2.5, shape = 22, fill = "white", color = y2c) +
  # geom_point(aes(y = y3), size = 2.5, shape = 24, fill = "white", color = y3c) +
  labs(x = "Time steps", y = "Prey2 to prey1 relative mean densities")

ggp <- fig + facet_grid(data$pry1cons~data$prdCtchProb1) +
  labs(x = "Time steps", y = "Prey2 to prey1 relative mean densities")

# save plot in this folder
ggsave(filename = paste("relativeDensities/relativeDensity-", XparamName, "X", YparamName, ".pdf", sep = ""), path = path, plot = ggp, width = 6.22, height = 5.73, limitsize = TRUE)
