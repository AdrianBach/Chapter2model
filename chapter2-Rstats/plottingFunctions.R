library(ggplot2)

path = "/home/adrian/Documents/GitKraken/Chapter2model/findEquilibrium/findEq-size25-res1max50-pry1init250-pry1Cons1-prdInit25-prdSurv10-prdCtch0.05/stats-findEq-size25-res1max50-pry1init250-pry1Cons1-prdInit25-prdSurv10-prdCtch0.05/statsResults-findEq-size25-res1max50-pry1init250-pry1Cons1-prdInit25-prdSurv10-prdCtch0.05.csv"

data <- read.csv(path)

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
fig

ggsave(filename = "test.pdf", plot = fig, path = "/home/adrian/Documents/GitKraken/Chapter2model/findEquilibrium/findEq-size25-res1max50-pry1init250-pry1Cons1-prdInit25-prdSurv10-prdCtch0.05/stats-findEq-size25-res1max50-pry1init250-pry1Cons1-prdInit25-prdSurv10-prdCtch0.05/", limitsize = TRUE)
