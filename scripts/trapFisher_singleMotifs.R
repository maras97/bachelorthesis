

## Read raw TRAP results of strong mESC Enhancers
mESC_high <- read.table("../STARRseq_enhancers/trapRaw_strongESC_singleMotifs.txt")
## Read raw TRAP results of weak mESC Enhancers
mESC_low <- read.table("../STARRseq_enhancers/trapRaw_weakESC_singleMotifs.txt")

## Read raw TRAP results of strong RARa Enhancers
RARa_high <- read.table("../STARRseq_enhancers/trapRaw_strongRARa_singleMotifs.txt")
## Read raw TRAP results of weak RARa Enhancers
RARa_low <- read.table("../STARRseq_enhancers/trapRaw_weakRARa_singleMotifs.txt")


## prepare value matrix for strong mESC enhancers
l_mESC <- nrow(mESC_high)/606  # 257
profiles <- mESC_high[1:606,4] 
values_mESC_high <- c()
for (k in 1:l_mESC){
  col <- mESC_high[(1+(k-1)*606):(k*606),5]
  values_mESC_high <- cbind(values_mESC_high,col)
}
rownames(values_mESC_high) <- profiles
colnames(values_mESC_high) <- seq(1,l_mESC,1)

## prepare value matrix for weak mESC enhancers
values_mESC_low <- c()
for (k in 1:l_mESC){
  col <- mESC_low[(1+(k-1)*606):(k*606),5]
  values_mESC_low <- cbind(values_mESC_low,col)
}
rownames(values_mESC_low) <- profiles
colnames(values_mESC_low) <- seq(1,l_mESC,1)

## prepare value matrix for strong RARa enhancers
l_RARa <- nrow(RARa_high)/606  # 1141
values_RARa_high <- c()
for (k in 1:l_RARa){
  col <- RARa_high[(1+(k-1)*606):(k*606),5]
  values_RARa_high <- cbind(values_RARa_high,col)
}
rownames(values_RARa_high) <- profiles
colnames(values_RARa_high) <- seq(1,l_RARa,1)

## prepare value matrix for weak RARa enhancers
values_RARa_low <- c()
for (k in 1:l_RARa){
  col <- RARa_low[(1+(k-1)*606):(k*606),5]
  values_RARa_low <- cbind(values_RARa_low,col)
}
rownames(values_RARa_low) <- profiles
colnames(values_RARa_low) <- seq(1,l_RARa,1)

## mESC
##
## APPLY A CUT-OFF OF K=200 TO SEPERATE TOP FROM RESIDUAL SEQUENCES 
## CREATE CONTINGENCY TABLE WITH VARIABLES TOP/RESIDUAL  VS. STRONG/WEAK 
## FISHER TESTS OF MESC ENHANCERS FOR EVERY MOTIF/PROFILE 
## ONE-SIDED FISHER-TEST WITH ALTERNATIVE = GREATER -> HIGH
k <- 200
pvalues_mESC_high <- c()
pvalues_mESC_low <- c()
for (p in 1:length(profiles)){
  tf <- append(values_mESC_high[p,],values_mESC_low[p,])
  tf <- cbind(tf,append(rep("high",l_mESC),rep("low",l_mESC)))
  tf <- tf[order(as.numeric(tf[,1]), decreasing = TRUE),]
  top_mESC_high <- sum(tf[1:k,2]=="high",na.rm=TRUE)
  cont <- rbind(c(top_mESC_high,k-top_mESC_high),c(l_mESC-top_mESC_high,l_mESC-(k-top_mESC_high)))
  pvalues_mESC_high <- append(pvalues_mESC_high,fisher.test(cont, alternative = "greater")$p.value)
  pvalues_mESC_low <- append(pvalues_mESC_low,fisher.test(cont, alternative = "less")$p.value)
}
## save results for strong mESC enhancers in txt-File
names(pvalues_mESC_high) <- profiles
pvalues_mESC_high <- sort(pvalues_mESC_high)
padj_mESC_high <- p.adjust(pvalues_mESC_high, method = "fdr")
write.table(cbind(pvalues_mESC_high,padj_mESC_high),"../RESULTS_trap/trap_strongESC_singleMotifs.txt", quote = F)
## save results for weak mESC enhancers in txt-File
names(pvalues_mESC_low) <- profiles
pvalues_mESC_low <- sort(pvalues_mESC_low)
padj_mESC_low <- p.adjust(pvalues_mESC_low, method = "fdr")
write.table(cbind(pvalues_mESC_low,padj_mESC_low),"../RESULTS_trap/trap_weakESC_singleMotifs.txt", quote = F)


## RARa
##
## APPLY A CUT-OFF OF K=200 TO SEPERATE TOP FROM RESIDUAL SEQUENCES 
## CREATE CONTINGENCY TABLE WITH VARIABLES TOP/RESIDUAL  VS. STRONG/WEAK 
## FISHER TESTS OF RARa ENHANCERS FOR EVERY MOTIF/PROFILE 
## ONE-SIDED FISHER-TEST WITH ALTERNATIVE = GREATER -> HIGH
k <- 500
pvalues_RARa_high <- c()
pvalues_RARa_low <- c()
for (p in 1:length(profiles)){
  tf <- append(values_RARa_high[p,],values_RARa_low[p,])
  tf <- cbind(tf,append(rep("high",l_RARa),rep("low",l_RARa)))
  tf <- tf[order(as.numeric(tf[,1]), decreasing = TRUE),]
  top_RARa_high <- sum(tf[1:k,2]=="high",na.rm=TRUE)
  cont <- rbind(c(top_RARa_high,k-top_RARa_high),c(l_RARa-top_RARa_high,l_RARa-(k-top_RARa_high)))
  pvalues_RARa_high <- append(pvalues_RARa_high,fisher.test(cont, alternative = "greater")$p.value)
  pvalues_RARa_low <- append(pvalues_RARa_low,fisher.test(cont, alternative = "less")$p.value)
}
## save results for strong RARa enhancers in txt-File
names(pvalues_RARa_high) <- profiles
pvalues_RARa_high <- sort(pvalues_RARa_high)
padj_RARa_high <- p.adjust(pvalues_RARa_high, method = "fdr")
write.table(cbind(pvalues_RARa_high,padj_RARa_high),"../RESULTS_trap/trap_strongRARa_singleMotifs.txt", quote = F)
## save results for weak RARa enhancers in txt-File
names(pvalues_RARa_low) <- profiles
pvalues_RARa_low <- sort(pvalues_RARa_low)
padj_RARa_low <- p.adjust(pvalues_RARa_low, method = "fdr")
write.table(cbind(pvalues_RARa_low,padj_RARa_low),"../RESULTS_trap/trap_weakRARa_singleMotifs.txt", quote = F)

