library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)


cl_id <- read.table("../motifDB/Clusters_1TF.tab", col.names = c("motif", "matrix", "TF"))


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
## 
plotHeatmapDiff <- function(mat) {
  mat <- mat[order(mat[,3],decreasing=T),]
  tf_id <- match(rownames(mat),cl_id[,1])
  names_cluster <- cl_id[tf_id,3]
  rownames(mat) <- paste(rownames(mat),names_cluster, sep = ":")
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
  tf_id <- match(rownames(mat),cl_id[,1])
  names_cluster <- cl_id[tf_id,3]
  rownames(mat) <- paste(rownames(mat),names_cluster, sep = ":")
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
motifcounter_esc_high <- read.table("../RESULTS_motifcounter/motifcounter_strongESC_clusteredMotifs.txt", fill = T, header = T, row.names = 1)
motifcounter_esc_low <- read.table("../RESULTS_motifcounter/motifcounter_weakESC_clusteredMotifs.txt", fill = T, header = T, row.names = 1)

# read motifcounter results for RARa 
motifcounter_rar_high <- read.table("../RESULTS_motifcounter/motifcounter_strongRARa_clusteredMotifs.txt", fill = T, header = T,row.names = 1)
motifcounter_rar_low <- read.table("../RESULTS_motifcounter/motifcounter_weakRARa_clusteredMotifs.txt", fill = T, header = T, row.names = 1)

## order every table according to alphabetical order of the TFBSs
order <- sort(rownames(motifcounter_esc_high))
motifcounter_esc_high <- motifcounter_esc_high[order,]
motifcounter_esc_low  <- motifcounter_esc_low[order,]
motifcounter_rar_high <- motifcounter_rar_high[order,]
motifcounter_rar_low  <- motifcounter_rar_low[order,]


## top ratios weak vs. strong ESC enhancers of foldChanges from Motifcounter
mat <- cbind(motifcounter_esc_low[,1],motifcounter_esc_high[,1])
mat <- cbind(mat,-log(mat[,1]/mat[,2]))
rownames(mat) <- order
## top ratios strong vs. weak
mat_top <- topSignif(mat, k=25)
colnames(mat_top) <- c("mESC\n weak","mESC\n strong","-log(Ratio)")
p <- plotHeatmapDiff(mat_top)
pdf("../heatmaps/motifcounter_mESC_strong_weak_topRatios_clusteredMotifs.pdf")
add_heatmap(p[[1]],p[[2]])
dev.off()

## top ratios RAR-weak vs. RAR-strong enhancer FoldChanges from Motifcounter
mat <- cbind(motifcounter_rar_low[,1],motifcounter_rar_high[,1])
mat <- cbind(mat,-log(mat[,1]/mat[,2]))
rownames(mat) <- order
## top ratios strong vs. weak
mat_top <- topSignif(mat, k=25)
colnames(mat_top) <- c("RARa\n weak","RARa\n strong","-log(Ratio)")
p <- plotHeatmapDiff(mat_top)
pdf("../heatmaps/motifcounter_RARa_strong_weak_topRatios_clusteredMotifs.pdf")
add_heatmap(p[[1]],p[[2]])
dev.off()

