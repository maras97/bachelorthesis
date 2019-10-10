library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)


## removes motifs without foldChange > 1 in one or both enhancer sets 
## and returns the 25 motifs with highest absolute values regarding log-ratio of foldChanges from motifcounter results
topSignif <- function(mat, k){
  signif_high <- which(mat[,1]>1) 
  signif_low <- which(mat[,2]>1) 
  signif <- unique(c(signif_high,signif_low)) 
  mat <- mat[signif,]
  top25 <- head(order(abs(mat[,3]), decreasing = T), n = k)
  mat <- mat[top25,]
  return(mat)
}
## creates a Heatmap with 2 columns for foldChanges from 1) weak enhancers 2) strong enhancers and a 3rd annotation column with log-ratio of FoldChanges
plotHeatmapDiff <- function(mat) {
  mat <- mat[order(mat[,3],decreasing=T),]
  pal <- c("indianred3", "white", "cadetblue")
  tw <- mat[,3]
  col_fun = colorRamp2(c(-max(abs(tw)), 0, max(abs(tw))),col = pal)
  ha = HeatmapAnnotation(logRatio =tw, name = "-log(Ratio)", 
                         col = list(logRatio = col_fun), which = "row",
                         width = unit(2.1, "cm"))
  h <- Heatmap(mat[,c(1,2)], cluster_rows = F, cluster_columns = F, row_title_side = "left", row_names_side = "left",
               row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), 
               col = colorRampPalette((brewer.pal(n = 9, name = "Purples")))(1000),
               column_names_max_height = unit(2.2, "cm"), name = "FoldChange")
  return(c(h,ha))
}
## creates a Heatmap of the input matrix with a given color scale without annotation column, only simple legend 
plotHeatmapScale <- function(mat) {
  tw <- c(mat)
  pal <- c("indianred3", "white", "cadetblue")
  col_fun = colorRamp2(c(-max(abs(tw)), 0, max(abs(tw))),col = pal)
  h <- Heatmap(mat, cluster_rows = F, cluster_columns = F, row_title_side = "left", row_names_side = "left",
               row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 12), 
               col = col_fun,
               column_names_max_height = unit(2.2, "cm"), name = "logRatio")
  return(h)
}
plotHeatmapScaleTRAP <- function(mat) {
  tw <- c(mat)
  pal <- c("indianred3", "white", "cadetblue")
  col_fun = colorRamp2(c(-max(abs(tw)), 0, max(abs(tw))),col = pal)
  h <- Heatmap(mat, cluster_rows = F, cluster_columns = F, row_title_side = "left", row_names_side = "left",
               row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 12), 
               col = col_fun,
               column_names_max_height = unit(2.2, "cm"), name = "logRatio")
  return(h)
}



# read motifcounter results for mESC 
motifcounter_esc_high <- read.table("../RESULTS_motifcounter/motifcounter_strongESC_singleMotifs.txt", fill = T, header = T, row.names = 1)
motifcounter_esc_low <- read.table("../RESULTS_motifcounter/motifcounter_weakESC_singleMotifs.txt", fill = T, header = T, row.names = 1)

# read motifcounter results for RARa 
motifcounter_rar_high <- read.table("../RESULTS_motifcounter/motifcounter_strongRARa_singleMotifs.txt", fill = T, header = T,row.names = 1)
motifcounter_rar_low <- read.table("../RESULTS_motifcounter/motifcounter_weakRARa_singleMotifs.txt", fill = T, header = T, row.names = 1)
# read TRAP results for RARa 
trap_rar_high <- read.table("../RESULTS_trap/trap_strongRARa_singleMotifs.txt", fill = T, header = T)
trap_rar_low <- read.table("../RESULTS_trap/trap_weakRARa_singleMotifs.txt", fill = T, header = T)


## order every table according to alphabetical order of the TFBSs
order <- sort(rownames(motifcounter_esc_high))
motifcounter_esc_high <- motifcounter_esc_high[order,]
motifcounter_esc_low  <- motifcounter_esc_low[order,]
motifcounter_rar_high <- motifcounter_rar_high[order,]
motifcounter_rar_low  <- motifcounter_rar_low[order,]

orderTRAP <- sort(rownames(trap_rar_high))
trap_rar_high <- trap_rar_high[orderTRAP,]
trap_rar_low  <- trap_rar_low[orderTRAP,]



## top ratios weak vs. strong ESC enhancers of foldChanges from Motifcounter
mat <- cbind(motifcounter_esc_low[,1],motifcounter_esc_high[,1])
mat <- cbind(mat,-log(mat[,1]/mat[,2]))
rownames(mat) <- order
mat_top <- topSignif(mat, k=25)
colnames(mat_top) <- c("mESC\n weak","mESC\n strong","-log(Ratio)")
p <- plotHeatmapDiff(mat_top)
pdf("../heatmaps/motifcounter_mESC_strong_weak_topRatios_singleMotifs.pdf")
add_heatmap(p[[1]],p[[2]])
dev.off()


## top ratios RAR-weak vs. RAR-strong enhancer FoldChanges from Motifcounter
mat <- cbind(motifcounter_rar_low[,1],motifcounter_rar_high[,1])
mat <- cbind(mat,-log(mat[,1]/mat[,2]))
rownames(mat) <- order
mat_top <- topSignif(mat, k=25)
colnames(mat_top) <- c("RARa\n weak","RARa\n strong","-log(Ratio)")
p <- plotHeatmapDiff(mat_top)
pdf("../heatmaps/motifcounter_RARa_strong_weak_topRatios_singleMotifs.pdf")
add_heatmap(p[[1]],p[[2]])
dev.off()


## comparison of RAR-RXR motifs for WEAK- vs. STRONG RARa enhancer foldChanges from Motifcounter
rar <- c("RARA-RXRA-DR0","RARA-RXRA-DR1","RARA-RXRA-DR2","RARA-RXRA-DR3","RARA-RXRA-DR4","RARA-RXRA-DR5","RARA-RXRA-DR6",
         "RARA-RXRA-DR7","RARA-RXRA-DR8",
         "RARA-RXRA-IR0","RARA-RXRA-IR1","RARA-RXRA-IR2","RARA-RXRA-IR3","RARA-RXRA-IR4","RARA-RXRA-IR5","RARA-RXRA-IR6",
         "RARA-RXRA-IR7","RARA-RXRA-IR8",
         "RARA-RXRA-ER0","RARA-RXRA-ER1","RARA-RXRA-ER2","RARA-RXRA-ER3","RARA-RXRA-ER4","RARA-RXRA-ER5","RARA-RXRA-ER6",
         "RARA-RXRA-ER7","RARA-RXRA-ER8", "RARA", "RARA(var.2)", "RARA-RXRG", "Rarb", "Rarb(var.2)", "Rarg","Rarg(var.2)")
mat <- cbind(motifcounter_rar_low[,1],motifcounter_rar_high[,1])
mat <- cbind(mat,-log(mat[,1]/mat[,2]))
rownames(mat) <- order
mat_rar <- mat[rar,]
colnames(mat_rar) <- c("RARa\n weak","RARa\n strong", "-log(Ratio)")
mat_rar <- mat_rar[order(abs(mat_rar[,3]), decreasing = T),]
p <- plotHeatmapDiff(mat_rar)
pdf("../heatmaps/motifcounter_RARa_strong_weak_RAmotifs.pdf")
add_heatmap(p[[1]],p[[2]])
dev.off()


## heatmap for comparison of RAR/RXR Motifs with varying spaces and orientation 
## with log-ratio of foldChanges from MOTIFCOUNTER between strong and weak RARa enhancers
mat_rar <- cbind(-log(motifcounter_rar_low[dr,1]/motifcounter_rar_high[dr,1]),
                 -log(motifcounter_rar_low[ir,1]/motifcounter_rar_high[ir,1]),
                 -log(motifcounter_rar_low[er,1]/motifcounter_rar_high[er,1]))
colnames(mat_rar) <- c("DR","IR", "ER")
rownames(mat_rar) <- seq(0,8)
pdf("../heatmaps/motifcounter_RARa_differential_3x8_RaraRxra.pdf")
plotHeatmapScale(mat_rar)
dev.off()


## heatmap for comparison of RAR/RXR Motifs with varying spaces and orientation 
## with log-ratio of values from TRAP between strong and weak RARa enhancers
mat_rar <- cbind(log(trap_rar_low[dr,1]/trap_rar_high[dr,1]),
                 log(trap_rar_low[ir,1]/trap_rar_high[ir,1]),
                 log(trap_rar_low[er,1]/trap_rar_high[er,1]))
colnames(mat_rar) <- c("DR","IR", "ER")
rownames(mat_rar) <- seq(0,8)
#write.table(signif(mat_rar), "motifcounter_RAR-RXR_Ratio.txt", quote = F, row.names = T, col.names = T)
pdf("../heatmaps/trap_RARa_differential_3x8_RaraRxra.pdf")
plotHeatmapScaleTRAP(mat_rar)
dev.off()
