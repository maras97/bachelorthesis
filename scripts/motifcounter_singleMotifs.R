library(motifcounter)
library(Biostrings)

## PREPROCESSING
##
## Read sequences of strong mESC Enhancers
mESC_file_high = "../STARRseq_enhancer/mESC_strongEnhancers.fasta"
mESC_high = readDNAStringSet(mESC_file_high)
## Read sequences of weak mESC Enhancers
mESC_low = readDNAStringSet(mESC_file_low)

## Read sequences of strong RARa Enhancers
RARa_file_high = "../STARRseq_enhancer/RARa_strongEnhancers.fasta"
RARa_high = readDNAStringSet(RARa_file_high)
## Read sequences of weak RARa Enhancers
RARa_file_low = "../STARRseq_enhancer/RARa_weakEnhancers.fasta"
RARa_low = readDNAStringSet(RARa_file_low)

## Union of both strong and weak enhancers for Background 
mESC_union = c(mESC_high,mESC_low)
RARa_union = c(RARa_high,RARa_low)

## Estimate order-1 background models
order = 1
bg_mESC = readBackground(mESC_union, order)
bg_RARa = readBackground(RARa_union, order)

## Load file with motifs in PSCM format 
motiffile = "../motifDB/jaspar_singleMotifs.pscm"
motifs <- read.table(motiffile, fill = T, stringsAsFactors = F)

## convert into list of motifs in PFM-format that can be used from the motifcounter package
names <- c()
M <- list()
c <- 1
# Iterate over Motifs in Motiffile to extract the matrices
for (m in 1:nrow(motifs)){
  if (substr(motifs[,1][m],1,1)==">"){
    tf <- gsub("/name=","",motifs[m,][2])
    names <- append(names,tf)
    p <- 1
    while((substr(motifs[,1][m+p],1,1)!=">") & ((m+p)<=nrow(motifs))){
      p <- p + 1
    }
    motif <- motifs[(m+1):(m+p-1),1:4]
    for (r in 1:nrow(motif)){
      motif[r,] = prop.table(data.matrix(motif[r,]))
    }
    motif <- t(data.matrix(motif))
    M[[c]] <- normalizeMotif(motif)
    c = c + 1
  }
}
names(M) <- names


## MOTIF ENRICHMENT ANALYIS mESC
##
## motifcounter analysis for strong mESC enhancer sequences
fold_mESC_high <- c()
pvalues_mESC_high <- c()
hits_mESC_high <- c()
## iteration over all motif profiles, motifcounter::motifEnrichment returns object with motif Enrichment information 
for (i in 1:length(M)){
  result_mESC_high = motifEnrichment(mESC_high, M[[i]], bg_mESC, singlestranded = FALSE, method = "compound")
  fold_mESC_high  <- append(fold_mESC_high,result_mESC_high$fold)
  pvalues_mESC_high  <- append(pvalues_mESC_high,result_mESC_high$pvalue)
  hits_mESC_high <- append(hits_mESC_high, sum(motifcounter:::numMotifHits(mESC_high, M[[i]], bg_mESC, singlestranded = FALSE)$numofhits))
}

# summarize resuls for strong mESC enhancers and save in txt-File 
names(pvalues_mESC_high) <- names(M)
padj_mESC_high <- p.adjust(pvalues_mESC_high, method = "fdr")
res_mESC_high <- cbind(names(M),fold_mESC_high,pvalues_mESC_high,padj_mESC_high,hits_mESC_high)
colnames(res_mESC_high) <- c("TFBS-name","fold-change","pvalue","adj-pvalue","number-of-hits")
res_mESC_high <- res_mESC_high[order(as.numeric(res_mESC_high[,3]), decreasing = F),]
write.table(res_mESC_high,"../RESULTS_motifcounter/motifcounter_strongESC_singleMotifs.txt", quote = F,row.names = F)


## motifcounter analysis for weak mESC enhancer sequences
fold_mESC_low <- c()
pvalues_mESC_low <- c()
hits_mESC_low <- c()
## iteration over all motif profiles, motifcounter::motifEnrichment returns object with motif Enrichment information 
for (i in 1:length(M)){
  result_mESC_low = motifEnrichment(mESC_low, M[[i]], bg_mESC, singlestranded = FALSE, method = "compound")
  fold_mESC_low  <- append(fold_mESC_low,result_mESC_low$fold)
  pvalues_mESC_low  <- append(pvalues_mESC_low,result_mESC_low$pvalue)
  hits_mESC_low <- append(hits_mESC_low, sum(motifcounter:::numMotifHits(mESC_low, M[[i]], bg_mESC, singlestranded = FALSE)$numofhits))
}

# summarize resuls for weak mESC enhancers and save in txt-File 
names(pvalues_mESC_low) <- names(M)
padj_mESC_low <- p.adjust(pvalues_mESC_low, method = "fdr")
res_mESC_low <- cbind(names(M),fold_mESC_low,pvalues_mESC_low,padj_mESC_low,hits_mESC_low)
colnames(res_mESC_low) <- c("TFBS-name","fold-change","pvalue","adj-pvalue","number-of-hits")
res_mESC_low <- res_mESC_low[order(as.numeric(res_mESC_low[,3]), decreasing = F),]
write.table(res_mESC_low,"../RESULTS_motifcounter/motifcounter_weakESC_singleMotifs.txt", quote = F,row.names = F)


## MOTIF ENRICHMENT ANALYIS RARa
##
##
## motifcounter analysis for strong RARa enhancer sequences
fold_RARa_high <- c()
pvalues_RARa_high <- c()
hits_RARa_high <- c()
## iteration over all motif profiles, motifcounter::motifEnrichment returns object with motif Enrichment information 
for (i in 1:length(M)){
  result_RARa_high = motifEnrichment(RARa_high, M[[i]], bg_RARa, singlestranded = FALSE, method = "compound")
  fold_RARa_high  <- append(fold_RARa_high,result_RARa_high$fold)
  pvalues_RARa_high  <- append(pvalues_RARa_high,result_RARa_high$pvalue)
  hits_RARa_high <- append(hits_RARa_high, sum(motifcounter:::numMotifHits(RARa_high, M[[i]], bg_RARa, singlestranded = FALSE)$numofhits))
}

# summarize resuls for strong RARa enhancers and save in txt-File 
names(pvalues_RARa_high) <- names(M)
padj_RARa_high <- p.adjust(pvalues_RARa_high, method = "fdr")
res_RARa_high <- cbind(names(M),fold_RARa_high,pvalues_RARa_high,padj_RARa_high,hits_RARa_high)
colnames(res_RARa_high) <- c("TFBS-name","fold-change","pvalue","adj-pvalue","number-of-hits")
res_RARa_high <- res_RARa_high[order(as.numeric(res_RARa_high[,3]), decreasing = F),]
write.table(res_RARa_high,"../RESULTS_motifcounter/motifcounter_strongRARa_singleMotifs.txt", quote = F,row.names = F)


## motifcounter analysis for weak RARa enhancer sequences
fold_RARa_low <- c()
pvalues_RARa_low <- c()
hits_RARa_low <- c()
## iteration over all motif profiles, motifcounter::motifEnrichment returns object with motif Enrichment information 
for (i in 1:length(M)){
  result_RARa_low = motifEnrichment(RARa_low, M[[i]], bg_RARa, singlestranded = FALSE, method = "compound")
  fold_RARa_low  <- append(fold_RARa_low,result_RARa_low$fold)
  pvalues_RARa_low  <- append(pvalues_RARa_low,result_RARa_low$pvalue)
  hits_RARa_low <- append(hits_RARa_low, sum(motifcounter:::numMotifHits(RARa_low, M[[i]], bg_RARa, singlestranded = FALSE)$numofhits))
}

# summarize resuls for weak RARa enhancers and save in txt-File 
names(pvalues_RARa_low) <- names(M)
padj_RARa_low <- p.adjust(pvalues_RARa_low, method = "fdr")
res_RARa_low <- cbind(names(M),fold_RARa_low,pvalues_RARa_low,padj_RARa_low,hits_RARa_low)
colnames(res_RARa_low) <- c("TFBS-name","fold-change","pvalue","adj-pvalue","number-of-hits")
res_RARa_low <- res_RARa_low[order(as.numeric(res_RARa_low[,3]), decreasing = F),]
write.table(res_RARa_low,"../RESULTS_motifcounter/motifcounter_weakRARa_singleMotifs.txt", quote = F,row.names = F)



