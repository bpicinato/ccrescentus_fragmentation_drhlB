######################
# bpicinato Apr 2023 #
######################

######################
#   DESCRIPTION     #
######################
# First part: script to calculate values of P for whole genome and generate plot
# Second part: calculate BER

######################

library("tidyverse")
library("ggplot2")

options(scipen=10000)

# sample 1 used for plot

# ATTENTION: using "threeprime" files that are actually fiveprime of transcript
# because of sequencing kit and protocol - dUTP protocol
# (Illumina Stranded Total RNA Prep Ligation with Ribozero Plus)
# also, strandness is inverted

##### iscR #####

# reading

##### inverted because of kit - 3 prime files, but actually real 5 prime #####

NA_fwd_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/iscR/sample1/threeprime/NA_all-3primeprofile-fwd.bedgraph", col_names = F)
NA_rev_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/iscR/sample1/threeprime/NA_all-3primeprofile-rev.bedgraph", col_names = F)
iscR_fwd_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/iscR/sample1/threeprime/iscR_all-3primeprofile-fwd.bedgraph", col_names = F)
iscR_rev_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/iscR/sample1/threeprime/iscR_all-3primeprofile-rev.bedgraph", col_names = F)

NA_fwd_3 = rep(c(NA_fwd_bedg_3$X4), times = NA_fwd_bedg_3$X3 - NA_fwd_bedg_3$X2)
NA_rev_3 = rep(c(NA_rev_bedg_3$X4), times = NA_rev_bedg_3$X3 - NA_rev_bedg_3$X2)

iscR_fwd_3 = rep(c(iscR_fwd_bedg_3$X4), times = iscR_fwd_bedg_3$X3 - iscR_fwd_bedg_3$X2)
iscR_rev_3 = rep(c(iscR_rev_bedg_3$X4), times = iscR_rev_bedg_3$X3 - iscR_rev_bedg_3$X2)

p_fwd_3_iscR = iscR_fwd_3 / (iscR_fwd_3 + NA_fwd_3)
p_rev_3_iscR = iscR_rev_3 / (iscR_rev_3 + NA_rev_3)
p_threeprime_iscR = c(p_fwd_3_iscR, p_rev_3_iscR)

#######

# putting into dataframe to plot
df = data.frame("P" = p_threeprime_iscR)
ggplot(df, aes(x = P)) + 
  geom_histogram(bins = 20, color = "white") + 
  geom_hline(yintercept = sum(abs(p_threeprime_iscR - 0) < 1e-6, na.rm = T), linetype = "dashed", 
              color = "grey", size = 0.25) + 
  theme_minimal() +
  labs(x="P", y = "Frequency") +
  ylim(0, 400000) +  
  theme(panel.border = element_rect(linetype = "solid", color = "black", fill = NA),
        legend.position = "none",
        text = element_text(color = "black", size = 8),
        axis.text = element_text(color = "black"), 
        strip.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10)
        )
# saving
ggsave("degradation_rhlB/final_plots/5prime-profile-iscR.png", 
       dpi = 600, width = 3, height = 3, units = "in", bg = "white")

################################################################################

###################
# tkoide Ago 2023 #
###################

#  Statistical analysis  #
# Bayes Error Rate - BER #
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC517707/
# Vencio et al., 2004

# getting info from histograms
h = hist(p_threeprime_iscR)
l = length(h$counts)

# calculating standard devitation without the extremes (P=0 and P=1)
sigma = sd(h$counts[c(-1, -l)])

# making error distribution for P=0 (Y0) and P=1 (Y1) with mean = (counts P=0) 
# or = (counts P=1) and sd = sigma
X = seq(0, max(h$counts[1] + 3*sigma, h$counts[l] + 3*sigma), len=2000)
Y0 = dnorm(X, h$counts[1], sigma)
Y1 = dnorm(X, h$counts[l], sigma)
# plotting dist
plot(X, Y0); lines(X, Y1)
# calculating BER -> area beneath the intersection of the curves * normalization
# see Fig. 1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC517707/figure/F1/
BER = sum(pmin(Y0, Y1)) * (X[2]-X[1])

# BER iscR: 0.2544957

