library(minfi)
library(data.table)
library(impute)
library(lumi)

#args = commandArgs(trailingOnly = TRUE)
#cancer = args[1]
#file = "GSE52865" # 276129*57
#file = "GSE72874-GPL13534" # 237020*200
#file = "GSE38266" # 276129*42
#file = "GSE61441" # 157528*92
#file = "GSE75041" # 276129*66
#file = "GSE97466" # 263522*141

zz = gzfile(paste0(file, '_beta.gz'),'rt')
mat = read.csv(zz, header=T)
#mat = fread(paste0("unzip -p ", file, "_beta.zip"))
matrix = data.frame(mat)[, -1]
row.names(matrix) = mat[, 1]
head(matrix)
dim(matrix) # 276916 201

## missing value imputation
sum(is.na(matrix))
matrix = as.matrix(matrix)
#matrix_nona = as.matrix(matrix)
matrix_nona = impute.knn(matrix, k = 10, rowmax = 0.25)$data
sum(is.na(matrix_nona))
head(matrix_nona)

## BMIQ
mset = makeGenomicRatioSetFromMatrix(matrix_nona) #276851
# origin.cpg = data.frame(ESCA_mat[, 1])$V1
# missing.cpg = unique(origin.cpg[! origin.cpg %in% ESCA.cpg])
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
beta.norm = apply(beta, 2, function(x) BMIQ(x, probe_types, nfit = 10000, plots=F)$nbeta)
#21 #29 #40 #42
head(beta.norm)

# beta.norm = apply(beta, 2, as.numeric)
# row.names(beta.norm) = row.names(beta)
                                 
## transfer Beta value to M value
#beta.norm = beta
m.norm = beta2m(beta.norm)

## output
dim(m.norm)
head(m.norm)
write.table(m.norm, file = gzfile(paste0(file, "_m_new.gz")), sep = " ")


### cancer purity
library(InfiniumPurify)

cancer = "BLCA"
mat = fread(cmd = paste0("unzip -p ", cancer, "_beta.zip"))
matrix = data.frame(mat[, !1])
row.names(matrix) = data.frame(mat[, 1])$V1
label = data.matrix(fread(paste0("label_", cancer, ".txt"), header = FALSE))
colnames(matrix) = label

tumor.data = matrix[, label == 2]
normal.data = matrix[, label == 1]

purity <- getPurity(tumor.data = tumor.data, normal.data = NULL, tumor.type= cancer)
tumor.purified = InfiniumPurify(tumor.data = tumor.data, normal.data = normal.data, 
                                purity = purity)

#purity.after <- getPurity(tumor.data = tumor.purified, normal.data = NULL, tumor.type= cancer)

