

sink("log/1_readData.log.txt", split=TRUE)

tryCatch({

# GBM

cat("Reading GBM data..")  
mat.GBM = data.matrix(read.csv("../data/GBM_counts.csv", header=TRUE, sep=",", 
                            stringsAsFactors=FALSE, check.names=FALSE))
  
pheno.GBM = read.csv("../data/GBM_pheno.csv", header=TRUE, sep=",")
cat("done.\n")

# PRAD

cat("Reading PRAD data..")  
mat.PRAD = data.matrix(read.csv("../data/PRAD_counts.csv", header=TRUE, sep=",", 
                               stringsAsFactors=FALSE, check.names=FALSE))

pheno.PRAD = read.csv("../data/PRAD_pheno.csv", header=TRUE, sep=",")
cat("done.\n")

# RCC

cat("Reading RCC data..")  
mat.RCC = data.matrix(read.csv("../data/RCC_counts.csv", header=TRUE, sep=",", 
                               stringsAsFactors=FALSE, check.names=FALSE))

pheno.RCC = read.csv("../data/RCC_pheno.csv", header=TRUE, sep=",")
cat("done.\n")

# BLCA

cat("Reading BLCA data..")  
mat.BLCA = data.matrix(read.csv("../data/BLCA_counts.csv", header=TRUE, sep=",", 
                               stringsAsFactors=FALSE, check.names=FALSE))

pheno.BLCA = read.csv("../data/BLCA_pheno.csv", header=TRUE, sep=",")
cat("done.\n")

mat.list = list(GBM=mat.GBM,
                PRAD=mat.PRAD, 
                RCC=mat.RCC,
                BLCA=mat.BLCA)

pheno.list = list(GBM=pheno.GBM,
                  PRAD=pheno.PRAD,
                  RCC=pheno.RCC,
                  BLCA=pheno.BLCA)


save(mat.list, pheno.list, file="../data/data.rda")

}, error = function(e){
  print(e)
})

print(sessionInfo())

rm(list=ls())
gc()

sink()


