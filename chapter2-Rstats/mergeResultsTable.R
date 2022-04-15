mergeResults <- function(path, keyword = c("Results", "Snapshot")) {
  
  # set work directory
  setwd(dir = path)
  
  # get the directory content
  content <- dir()
  
  # order alphabetically
  content <- content[order(content)]
  
  # only the csv files
  content <- grep(pattern = c(".csv"), x = content, value = T)
  
  # select the right results type ("Results or Snapshot")
  content <- grep(pattern = keyword, x = content, value = T)
  
  # read one file to get column number and headers
  zz <- read.csv(content[1])
  
  # create an empty table 
  tab <- data.frame(matrix(ncol = dim(zz)[2] + 1, nrow = 0))
  colnames(tab) <- c("replicate", colnames(zz))
  rm(zz)  # remove file from ram
  
  # loop over the files
  for (i in 1:length(content)) {
    
    # read file
    res <- read.csv(content[i])
    
    # bind the replicate number to the table
    replicate <- rep(i, times = dim(res)[1])
    res <- cbind(replicate, res)
    
    # append to the merge table
    tab <- rbind(tab, res)
  }
  
  # return table
  return(tab)
}
