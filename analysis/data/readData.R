# read individual data files and create full data set for each experiment
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

setwd(dir = 'analysis/data')

# function to read in datasets
readData <- function(dir, len = 540){
  # function to read in data files
  # dir = path to files
  # len = how long should each file be (rows)
  files <- list.files(path = dir, pattern = 'Participant')
  tmp.df <- data.frame()
  for (f in files){
    tmp <- read.csv(file = paste0(dir, f))
    tmp$ID <- f
    if (nrow(tmp) == len){ # full data set
      tmp.df <- rbind(tmp.df, tmp)
    } else{
      print(paste(f, 'is not complete'))
    }
  }
  return(tmp.df)
}

# read in
exp1data <- readData(dir = 'exp1/')
exp2data <- readData(dir = 'exp2/', len = 544)
exp3data <- readData(dir = 'exp3/')
exp4data <- readData(dir = 'exp4/')

# lets look at performance across participants

exp1data$Acc <- with(exp1data, as.integer(testChange == respChange))
exp2data$Acc <- with(exp2data, as.integer(testChange == respChange))
exp3data$Acc <- with(exp3data, as.integer(testChange == respChange))
exp4data$Acc <- with(exp4data, as.integer(testChange == respChange))

# participants performing lower than 60% are excluded
exp1pptAcc <- aggregate(Acc~ID, data = exp1data, FUN = mean)
exp1Excl <- subset(exp1pptAcc, Acc < .6) # 1 participant

exp1data <- subset(exp1data, !(ID %in% exp1Excl$ID))

exp2pptAcc <- aggregate(Acc~ID, data = exp2data, FUN = mean)
subset(exp2pptAcc, Acc < .6) # none

exp3pptAcc <- aggregate(Acc~ID, data = exp3data, FUN = mean)
subset(exp3pptAcc, Acc < .6) # none

exp4pptAcc <- aggregate(Acc~ID, data = exp4data, FUN = mean)
exp4Excl <- subset(exp4pptAcc, Acc < .6) # 2

exp4data <- subset(exp4data, !(ID %in% exp4Excl$ID))

# create participant numbers

exp1data$ppt <- as.numeric(as.factor(exp1data$ID))
exp2data$ppt <- as.numeric(as.factor(exp2data$ID))
exp3data$ppt <- as.numeric(as.factor(exp3data$ID))
exp4data$ppt <- as.numeric(as.factor(exp4data$ID))

# write files

write.table(exp1data, file = 'exp1/exp1_full')
write.table(exp2data, file = 'exp2/exp2_full')
write.table(exp3data, file = 'exp3/exp3_full')
write.table(exp4data, file = 'exp4/exp4_full')

