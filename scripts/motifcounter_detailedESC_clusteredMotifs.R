library(Biostrings)
library(motifcounter)
library(foreach)
library(doParallel)
library(doSNOW)

##PREPROCESSING
##
## Read sequences of strong mESC Enhancers
mESC_file_high = "../STARRseq_enhancer/mESC_strongEnhancers.fasta"
mESC_high = readDNAStringSet(mESC_file_high)
## Read sequences of weak mESC Enhancers
mESC_file_low = "../STARRseq_enhancer/mESC_weakEnhancers.fasta"
mESC_low = readDNAStringSet(mESC_file_low)


## Union of both strong and weak enhancers for Background 
mESC_union = c(mESC_high,mESC_low)
## Estimate order-1 background models
order = 1
bg_mESC = readBackground(mESC_union, order)

## load RDS-file with previously created list of PFMs of clustered motifs
M <- readRDS("../motifDB/motifs_cluster.rds")

## DETAILED MOTIF ENRICHMENT ANALYSIS WITH MOTIFCOUTNER FOR STRONG mESC ENHANCERS 
## set parameters for parallel processing 
registerDoParallel(32)
cl <- makeCluster(32)
registerDoSNOW(cl)
## show progress bar while iterating over the tasks
iterations <- length(mESC_high)*length(M)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
## store motif enrichment result for every sequence and motif individually 
result_mESC_high <- foreach(i = 1:length(mESC_high), .combine='rbind', .packages='motifcounter',.options.snow = opts) %:% 
  foreach(j = 1:length(M) , .combine='c', .packages='motifcounter') %dopar% {
    motifEnrichment(mESC_high[i], M[[j]], bg_top, singlestranded = FALSE, method = "compound")$fold
  }
close(pb)
stopCluster(cl)

## save result in a RDS-file 
saveRDS(result_mESC_high,"../RESULTS_motifcounter/motifcounter_strongESC_detailed_clusteredMotifs.rds")
#result_mESC_high <- readRDS("RESULTS_motifcounter/motifcounter_strongESC_detailed_clusteredMotifs.rds")


## DETAILED MOTIF ENRICHMENT ANALYSIS WITH MOTIFCOUTNER FOR STRONG mESC ENHANCERS 
## set parameters for parallel processing 
cl <- makeCluster(32)
registerDoSNOW(cl)
## show progress bar while iterating over the tasks
iterations <- length(mESC_low)*length(M)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
## store motif enrichment result for every sequence and motif individually 
result_mESC_low <- foreach(i = 1:length(mESC_low), .combine='rbind', .packages='motifcounter',.options.snow = opts) %:% 
  foreach(j = 1:length(M) , .combine='c', .packages='motifcounter') %dopar% {
    motifEnrichment(mESC_low[i], M[[j]], bg_top, singlestranded = FALSE, method = "compound")$fold
  }
close(pb)
stopCluster(cl)

## save result in a RDS-file 
saveRDS(result_mESC_low,"../RESULTS_motifcounter/motifcounter_weakESC_detailed_clusteredMotifs.rds")
#result_mESC_low <- readRDS("RESULTS_motifcounter/motifcounter_weakESC_detailed_clusteredMotifs.rds")


## combine results for strong and weak mESC enhancers in a matrix with information about enhancer activity and enrichment of any motif cluster
outMatrix <- rbind(result_mESC_high,result_mESC_low)
rownames(outMatrix) <- c(sapply(r, function(x) paste("strong",x)),sapply(r, function(x) paste("weak",x)))
colnames(outMatrix) <- names(M)


## OUTPUT
## save combined matrix in a RDS-file and additionally in a txt-file
saveRDS(outMatrix,"../RESULTS_motifcounter/motifcounter_StrongWeakESC_detailed_clusteredMotifs.rds")
write.table(outMatrix, "../RESULTS_motifcounter/motifcounter_StrongWeak_detailed_clusteredMotifs.txt", quote = F, row.names = T, col.names = T)




