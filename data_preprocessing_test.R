library(minfi)
library(data.table)
library(impute)
library(lumi)

args = commandArgs(trailingOnly = TRUE)
file = args[1]
#file = "GSE52865" # 276129*57
#file = "GSE72874-GPL13534" # 237020*200
#file = "GSE38266" # 276129*42
#file = "GSE61441" # 157528*92
#file = "GSE75041" # 276129*66
#file = "GSE97466" # 263522*141

## unzip all gz
library(GEOquery)
idatFiles <- list.files(paste0("~/Desktop/",file,"_RAW/"), pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

## read index
probe_names = read.csv(paste0("~/Desktop/",file,"_RAW/all_probe_names.txt"), header = FALSE)$V1 #276129

## read idat data
rgSet <- read.metharray.exp(paste0("~/Desktop/",file,"_RAW/"))
rgSet
head(sampleNames(rgSet))
getManifest(rgSet)
MSet <- preprocessRaw(rgSet) 
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet
beta <- getBeta(RSet)

## selecct data
sel_beta = beta[probe_names, ] #276129*36

## missing value imputation
sum(is.na(sel_beta))
matrix_nona = impute.knn(sel_beta, k = 10, rowmax = 0.25)$data

## BMIQ
mset = makeGenomicRatioSetFromMatrix(matrix_nona) #276129
probe_types = getProbeType(mset)
for(i in 1:length(probe_types)){
  if(isTRUE(probe_types[i] == "I")){
    probe_types[i] = "1"
  } else if (isTRUE(probe_types[i] == "II")){
    probe_types[i] = "2"
  }
}

options(connectionObserver = NULL)
library(wateRmelon)
beta = getBeta(mset)
beta.norm = apply(beta, 2, function(x) BMIQ(x, probe_types, nfit = 10000)$nbeta)

## transfer Beta value to M value
m.norm = beta2m(beta.norm)

## output
write.table(m.norm, file = gzfile(paste0(file, "_m.gz")), sep = " ")
