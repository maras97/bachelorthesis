#!/bin/bash


## TRAP for high- and low-scored sequences from mESC-STARR-enhancers, regarding every SINGLE MOTIF
./TRAP -s /../STARRseq_enhancer/mESC_strongEnhancers.fasta \
       -matrix /../motifDB/jaspar_singleMotifs.pscm \
       -thread 20  > /../RESULTS_trap/trapRaw_strongESC_singleMotifs.txt

./TRAP -s /../STARRseq_enhancer/mESC_weakEnhancers.fasta \
        -matrix /../motifDB/jaspar_singleMotifs.pscm \
        -thread 20  > /../RESULTS_trap/trapRaw_weakESC_singleMotifs.txt

## TRAP for high- and low-scored sequences from RARa-STARR-enhancers, regarding every SINGLE MOTIF
./TRAP -s /../STARRseq_enhancer/RARa_strongEnhancers.fasta \
        -matrix /../motifDB/jaspar_singleMotifs.pscm \
        -thread 20  > /../RESULTS_trap/trapRaw_strongRARa_singleMotifs.txt

./TRAP -s /../STARRseq_enhancer/RARa_weakEnhancers.fasta \
        -matrix /../motifDB/jaspar_singleMotifs.pscm \
        -thread 20  > /../RESULTS_trap/trapRaw_weakRARa_singleMotifs.txt


## TRAP for high- and low-scored sequences from mESC-STARR-enhancers, regarding CLUSTERED MOTIFS
./TRAP -s /../STARRseq_enhancer/mESC_strongEnhancers.fasta \
       -matrix /../motifDB/jaspar_clusteredMotifs.pscm \
       -thread 20  > /../RESULTS_trap/trapRaw_strongESC_clusteredMotifs.txt

./TRAP -s /../STARRseq_enhancer/mESC_weakEnhancers.fasta \
        -matrix /../motifDB/jaspar_clusteredMotifs.pscm \
        -thread 20  > /../RESULTS_trap/trapRaw_weakESC_clusteredMotifs.txt

## TRAP for high- and low-scored sequences from RARa-STARR-enhancers, regarding CLUSTERED MOTIFS
./TRAP -s /../STARRseq_enhancer/RARa_strongEnhancers.fasta \
        -matrix /../motifDB/jaspar_clusteredMotifs.pscm \
        -thread 20  > /../RESULTS_trap/trapRaw_strongRARa_clusteredMotifs.txt

./TRAP -s /../STARRseq_enhancer/RARa_weakEnhancers.fasta \
        -matrix /../motifDB/jaspar_clusteredMotifs.pscm \
        -thread 20  > /../RESULTS_trap/trapRaw_weakRARa_clusteredMotifs.txt
