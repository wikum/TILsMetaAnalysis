


tryCatch({

dir.create("tables")
dir.create("fig")
dir.create("obj")

setwd("analysis/")

dir.create("log")

source("1_readData.R")

rm(list=ls())
gc()
  
source("2_analysis.R")

rm(list=ls())
gc()

source("3_analysis.R")

rm(list=ls())
gc()

source("4_analysis.R")
  
}, error = function(e){
  print(e)
})

print(sessionInfo())

rm(list=ls())
gc()


