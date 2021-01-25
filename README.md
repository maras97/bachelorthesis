this repository contains scripts and data used for motif enrichment analysis experimentally identified enhancers by starr-seq

overview of folders: 

scripts: 
    r-scripts and bash-scripts used for the analysis 
   
motifdb: 
    pscm-matrices from the jaspar vertebrates core collection and jaspar core matrix clustering;
    additionally a table with an overview of single-tf motifs associated with the motif clusters
    
starrsseq_enhancer:  
    bed-file with genomic regions of active mesc enhancers obtained by starrseq, filtered for promoters and width-unified;
    additionally fasta-files with genomic sequences of weak and strong mesc and rara enhancers
    
results_motifcounter: 
    tables with results of motif enrichment analysis by motifcounter on weak and strong mesc and rara enhancers with single and clustered motifs;
    in each table motifs are ordered in descending order by their enrichment;
    additionally detailed motifcounter results for mesc enhancers with clustered motifs in txt- and rds-file to use for classifier later on 
        
results_trap: 
    tables with results of motif enrichment analysis by trap on weak and strong mesc and rara enhancers with single and clustered motifs;
    in each table motifs are ordered in descending order by their enrichment (in ascending order by p-value)
    
classifier: 
    rds-files with test- and training-sets and glmnet-object to reproduce resutls if needed;
    additionally table with beta-coefficients of the trained regression model
