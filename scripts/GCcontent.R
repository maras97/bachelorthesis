library("Biostrings")
library("ggplot2")
library("ggpubr")


## mESC
##
# read fasta-file with mESCenhancer sequences with strongly active enhancers
high_seq_ESC = readDNAStringSet("../STARRseq_enhancer/mESC_strongEnhancers.fasta")
# read fasta-file with mESC enhancer sequences with weakly active enhancers
low_seq_ESC = readDNAStringSet("./STARRseq_enhancer/mESC_weakEnhancers.fasta")

# get unified mESC enhancer sequence width
l_ESC <- high_seq_ESC@ranges@width[1]

# calculate ratio of CG to total number of nucleotides (enhancer width l) for every sequence in either strong or weak mESC enhancer set
ratio_high_ESC <- letterFrequency(high_seq_ESC, letters="CG")/l_ESC
ratio_low_ESC <- letterFrequency(low_seq_ESC, letters="CG")/l_ESC

# print mean GC-content in both strong and weak mESC enhancers in percent 
print(paste("strong mESC Enhancers mean CG-content: ",signif(mean(ratio_high_ESC)*100),"%",sep = ""))
print(paste("weak mESC Enhancers mean CG-content:",signif(mean(ratio_low_ESC)*100),"%",sep = ""))



## prepare ESC data frames to use in ggplot function for boxplots and wilcoxon test
GC_content <- c(ratio_low_ESC, ratio_high_ESC)
enhancer <- c(rep("weak",257), rep("strong",257))
state <- rep("mESC",514)
df <- data.frame(GC_content, enhancer, state)


## RARa
##
# read fasta-file with RARa enhancer sequences with strongly active enhancers
high_seq_RA = readDNAStringSet("RARa_strongEnhancers.fasta")
# read fasta-file with RARa enhancer sequences with weakly active enhancers
low_seq_RA = readDNAStringSet("RARa_weakEnhancers.fasta")

# get unified RARa enhancer sequence width
l_RA <- high_seq_RA@ranges@width[1]

# calculate ratio of CG to total number of nucleotides (enhancer width l) for every sequence in either strong or weak mESC enhancer set
ratio_high_RA <- letterFrequency(high_seq_RA, letters="CG")/l_RA
ratio_low_RA <- letterFrequency(low_seq_RA, letters="CG")/l_RA

# print mean GC-content in both strong and weak RARa enhancers in percent 
print(paste("strong RARa Enhancers mean CG-content: ",signif(mean(ratio_high_RA)*100),"%",sep = ""))
print(paste("weak RARa Enhancers mean CG-content:",signif(mean(ratio_low_RA)*100),"%",sep = ""))

## prepare RA data frames to use in ggplot function for boxplots and wilcoxon test
GC_content <- c(ratio_low_RA, ratio_high_RA)
enhancer <- c(rep("weak",1141), rep("strong",1141))
state <- rep("RARa",2282)
df_RA <- data.frame(GC_content, enhancer, state)
df <- rbind(df, df_RA)


## produce boxplot for comparison of GC-content distribution in weak vs. strong enhancers of mESC and rARa enhancers
## using ggplot functions
pdf("../GCcontent_boxplot.pdf")
p <- ggboxplot(df, x = "enhancer", y = "GC_content",
               color = "enhancer", palette = c("#e0ce9b", "#918050"), ylab = "GC - content", xlab = "",
               facet.by = "state", short.panel.labs = T, order = c("weak", "strong"), lwd = 0.8)
## add Wilcoxon - Test p-values for comparison of distributions
p + stat_compare_means(method.args = list(alternative = "two.sided"),
                       label.x = 1.1,  label.y = 0.79)
dev.off()

