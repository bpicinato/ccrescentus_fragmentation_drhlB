# RNA fragmentation analysis in *Caulobacter crescentus* using RNA-Seq data

Repository for scripts used in assessing fragmentation excess in  *Caulobacter crescentus* NA1000 *rhlB* mutant at 10°C and 30°C.

First script ramdomly samples n reads from BAM alingment files and generate 5' and 3' profiles as in https://github.com/alanlorenzetti/frtc. This was used to guarantee all libraries had the same number of alinged reads for further analysis. Second script calculated the proportion P for all genomic positions using the undersampled files and plots the resulting histograms.

This repository is part of study conducted primarily by Hugo Libonati de Araújo and Dr. Marilis do Valle Marques.
