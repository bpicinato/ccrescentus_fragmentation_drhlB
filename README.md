# RNA fragmentation analysis in *Caulobacter crescentus* using RNA-Seq data

Repository for scripts used in assessing fragmentation excess in the article "The DEAD-box RNA helicase RhlB is required for efficient RNA processing at
low temperature in Caulobacter".

Scripts used in calculating fragmentation excess in *Caulobacter crescentus* NA1000 *rhlB* mutant at 10°C and 30°C (`scripts/rhlB` - Figure 6B); *C. crescentus* NA1000 *iscR* mutant (`scripts/iscR` - Figure S6); *Pseudomonas aeruginosa* mutants *rhlE1*, *rhlE2*, and *rhlE1/rhlE2* ([Hausmann *et al.*, 2021](https://doi.org/10.1093/nar/gkab503)) (`scripts/rhlE` - Figure 6C). Each folder inside `scripts/` include 01) script for undersampling analysis and 02) script for calculating P, plotting the histogram and calculating Bayes Error Rate (BER). Scripts for *rhlB* and *iscR* mutants are pratically the same, changing only samples names and conditions. Script 01 for *rhlE* mutants is different from the others because it was made for single-end sequencing data.

Briefly, `01-undersampling.sh` ramdomly samples n reads (n being the amount of reads uniquely aligned in the library with less aligned reads among the samples) from BAM alingment files and generate 5' and 3' profiles as in https://github.com/alanlorenzetti/frtc; this was used to guarantee all libraries had the same number of alinged reads for further analysis. Script `02-calculatingP_plotting.R` calculates the proportion P for all genomic positions using the undersampled files, plots the resulting histograms and calculates the BER.

This repository is part of study conducted primarily by Hugo Libonati de Araújo and Dr. Marilis do Valle Marques.

If you use any content of this repository, please cite [The DEAD-box RNA helicase RhlB is required for efficient RNA processing at low temperature in *Caulobacter*](https://journals.asm.org/doi/10.1128/spectrum.01934-23).
