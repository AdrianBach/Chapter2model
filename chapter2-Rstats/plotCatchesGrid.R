x = data$timeStep
y1 = data$prey1catchesMean
y2 = data$prey2catchesMean
# y3 = data$predator1PopulationSizeMean
y1min = data$prey1catchesICinf
y2min = data$prey2catchesICinf
# y3min = data$predator1PopulationSizeICinf
y1max = data$prey1catchesICsup
y2max = data$prey2catchesICsup
# y3max = data$predator1PopulationSizeICsup
y1c = "red"
y2c = "blue"
# y3c = "orange"
tIntro = 200

fig <- ggplot(data, aes(x)) + 
  geom_rect(aes(xmin = 0, xmax = tIntro, ymin = 0, ymax = 1.05*max(y1max)), alpha=0.5, fill = "lightgrey") +
  geom_ribbon(aes(ymin = y1min, ymax = y1max), alpha = 0.2, size = 0.1, col = y1c, fill = y1c) +
  geom_ribbon(aes(ymin = y2min, ymax = y2max), alpha = 0.2, size = 0.1, col = y2c, fill = y2c) +
  # geom_ribbon(aes(ymin = y3min, ymax = y3max), alpha = 0.2, size = 0.1, col = y3c, fill = y3c) +
  geom_line(aes(y = y1), color = y1c, size = 0.5) +
  geom_line(aes(y = y2), color = y2c, size = 0.5) # +
  # geom_line(aes(y = y3), color = y3c)

ggp <- fig + facet_grid(data$pry1cons~data$prdCtchProb1) +
  labs(x = "Time steps", y = "Nb of catches")

# save plot in this folder
ggsave(filename = paste("catches-", XparamName, "X", YparamName, ".pdf", sep = ""), path = dirPath, plot = ggp, width = 6.22, height = 5.73, limitsize = TRUE)
