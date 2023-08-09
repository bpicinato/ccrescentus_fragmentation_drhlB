######################
# bpicinato Jul 2023 #
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

# ATTENTION: using "threeprime" files that are actually fiveprime of transcript
# because of sequencing kit and protocol - dUTP protocol
# (Illumina TruSeq stranded mRNA kit)
# also, strandness is inverted

# reading

##### inverted because of kit - 3 prime files, but actually real 5 prime #####
## 4 samples
# WT
# drhlE1
# drhlE2
# drhlE12

# WT
WT_fwd_bedg_3 = read_tsv("degradation_rhlB/rhlE_Paer/input_undersampling/threeprime/WT-3primeprofile-fwd.bedgraph", col_names = F)
WT_rev_bedg_3 = read_tsv("degradation_rhlB/rhlE_Paer/input_undersampling/threeprime/WT-3primeprofile-rev.bedgraph", col_names = F)
WT_fwd_3 = rep(c(WT_fwd_bedg_3$X4), times = WT_fwd_bedg_3$X3 - WT_fwd_bedg_3$X2)
WT_rev_3 = rep(c(WT_rev_bedg_3$X4), times = WT_rev_bedg_3$X3 - WT_rev_bedg_3$X2)
#rhlE1
RhlE1_fwd_bedg_3 = read_tsv("degradation_rhlB/rhlE_Paer/input_undersampling/threeprime/RhlE1-3primeprofile-fwd.bedgraph", col_names = F)
RhlE1_rev_bedg_3 = read_tsv("degradation_rhlB/rhlE_Paer/input_undersampling/threeprime/RhlE1-3primeprofile-rev.bedgraph", col_names = F)
RhlE1_fwd_3 = rep(c(RhlE1_fwd_bedg_3$X4), times = RhlE1_fwd_bedg_3$X3 - RhlE1_fwd_bedg_3$X2)
RhlE1_rev_3 = rep(c(RhlE1_rev_bedg_3$X4), times = RhlE1_rev_bedg_3$X3 - RhlE1_rev_bedg_3$X2)
#rhlE2
RhlE2_fwd_bedg_3 = read_tsv("degradation_rhlB/rhlE_Paer/input_undersampling/threeprime/RhlE2-3primeprofile-fwd.bedgraph", col_names = F)
RhlE2_rev_bedg_3 = read_tsv("degradation_rhlB/rhlE_Paer/input_undersampling/threeprime/RhlE2-3primeprofile-rev.bedgraph", col_names = F)
RhlE2_fwd_3 = rep(c(RhlE2_fwd_bedg_3$X4), times = RhlE2_fwd_bedg_3$X3 - RhlE2_fwd_bedg_3$X2)
RhlE2_rev_3 = rep(c(RhlE2_rev_bedg_3$X4), times = RhlE2_rev_bedg_3$X3 - RhlE2_rev_bedg_3$X2)
#rhlE12
RhlE12_fwd_bedg_3 = read_tsv("degradation_rhlB/rhlE_Paer/input_undersampling/threeprime/RhlE12-3primeprofile-fwd.bedgraph", col_names = F)
RhlE12_rev_bedg_3 = read_tsv("degradation_rhlB/rhlE_Paer/input_undersampling/threeprime/RhlE12-3primeprofile-rev.bedgraph", col_names = F)
RhlE12_fwd_3 = rep(c(RhlE12_fwd_bedg_3$X4), times = RhlE12_fwd_bedg_3$X3 - RhlE12_fwd_bedg_3$X2)
RhlE12_rev_3 = rep(c(RhlE12_rev_bedg_3$X4), times = RhlE12_rev_bedg_3$X3 - RhlE12_rev_bedg_3$X2)


## 3 comparisons
# RhlE1
p_fwd_3_1 = RhlE1_fwd_3 / (RhlE1_fwd_3 + WT_fwd_3)
p_rev_3_1 = RhlE1_rev_3 / (RhlE1_rev_3 + WT_rev_3)
p_threeprime_1 = c(p_fwd_3_1, p_rev_3_1)
# RhlE2
p_fwd_3_2 = RhlE2_fwd_3 / (RhlE2_fwd_3 + WT_fwd_3)
p_rev_3_2 = RhlE2_rev_3 / (RhlE2_rev_3 + WT_rev_3)
p_threeprime_2 = c(p_fwd_3_2, p_rev_3_2)
# RhlE12
p_fwd_3_12 = RhlE12_fwd_3 / (RhlE12_fwd_3 + WT_fwd_3)
p_rev_3_12 = RhlE12_rev_3 / (RhlE12_rev_3 + WT_rev_3)
p_threeprime_12 = c(p_fwd_3_12, p_rev_3_12)


#######

# putting into dataframe to plot
df = data.frame("RhlE1" = p_threeprime_1, "RhlE2" = p_threeprime_2, "RhlE12" = p_threeprime_12)
df_pivot = df %>% 
  pivot_longer(everything(), names_to = "mutant", values_to = "P" )

# making dataframe with horizontal line information for plot
lines = data.frame(mutant = c("RhlE1", "RhlE2", "RhlE12"), 
                   line = c(sum(abs(p_threeprime_1 - 0) < 1e-6, na.rm = T), 
                            sum(abs(p_threeprime_2 - 0) < 1e-6, na.rm = T),
                            sum(abs(p_threeprime_12 - 0) < 1e-6, na.rm = T)))
# labels for facets 
# https://community.rstudio.com/t/subscripts-and-superscripts-facet-wrap-facet-labels/81072/6
labels = as_labeller(c(RhlE1 = "Delta*italic(rhlE1)^Pa", 
                       RhlE2 = "Delta*italic(rhlE2)^Pa",
                       RhlE12 = "Delta*italic(rhlE1)^Pa/Delta*italic(rhlE2)^Pa"), default = label_parsed)
# plotting
ggplot(df_pivot, aes(x = P, fill = factor(mutant))) + 
  geom_histogram(bins = 20, color = "white") + 
  facet_wrap(~mutant, labeller = labels) +
  geom_abline(data = lines, aes(intercept = line, slope = 0), linetype = "dashed", 
              color = "BLACK", size = 0.25) + 
  theme_minimal() +
  scale_fill_manual(values = c("#B5CEA1", "#81AB5F", "#59783F")) +
  labs(x="P", y = "Frequency") +
  theme(panel.border = element_rect(linetype = "solid", color = "black", fill = NA),
        legend.position = "none",
        text = element_text(color = "black", size = 8, family="sans"),
        axis.text = element_text(color = "black", size = 8),
        strip.text.x = element_text(size = 10, family="sans"))
# saving
ggsave("degradation_rhlB/rhlE_Paer/5prime-profile_rhlE.png", 
       dpi = 300, width = 15, height = 7.5, units = "cm", bg = "white")

################################################################################

###################
# tkoide Ago 2023 #
###################

#  Statistical analysis  #
# Bayes Error Rate - BER #
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC517707/
# Vencio et al., 2004

# getting info from histograms
h1 = hist(p_threeprime_1)
l1 = length(h1$counts)
h2 = hist(p_threeprime_2)
l2 = length(h2$counts)
h12 = hist(p_threeprime_12)
l12 = length(h12$counts)

# calculating standard devitation without the extremes (P=0 and P=1)
sigma1 = sd(h1$counts[c(-1, -l1)])
sigma2 = sd(h2$counts[c(-1, -l2)])
sigma12 = sd(h12$counts[c(-1, -l12)])

## RHlE1 ##
# making error distribution for P=0 (Y0) and P=1 (Y1) with mean = (counts P=0) 
# or = (counts P=1) and sd = sigma
X1 = seq(0, max(h1$counts[1] + 3*sigma1, h1$counts[l1] + 3*sigma1), len=2000)
Y0_1 = dnorm(X1, h1$counts[1], sigma1)
Y1_1 = dnorm(X1, h1$counts[l1], sigma1)
# plotting dist
plot(X1, Y0_1); lines(X1, Y1_1)
# calculating BER -> area beneath the intersection of the curves * normalization
# see Fig. 1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC517707/figure/F1/
BER = sum(pmin(Y0_1, Y1_1)) * (X1[2]-X1[1])

# BER rhlE1: 0.9171446

## RhlE2 ##
# making error distribution for P=0 (Y0) and P=1 (Y1) with mean = (counts P=0) 
# or = (counts P=1) and sd = sigma
X2 = seq(0, max(h2$counts[1] + 3*sigma2, h2$counts[l2] + 3*sigma2), len=2000)
Y0_2 = dnorm(X2, h2$counts[1], sigma2)
Y1_2 = dnorm(X2, h2$counts[l2], sigma2)
# plotting dist
plot(X2, Y0_2); lines(X2, Y1_2)
# calculating BER -> area beneath the intersection of the curves * normalization
# see Fig. 1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC517707/figure/F1/
BER = sum(pmin(Y0_2, Y1_2)) * (X2[2]-X2[1])

# BER rhlE2: 0.0003501258

## RHlE12 ##
# making error distribution for P=0 (Y0) and P=1 (Y1) with mean = (counts P=0) 
# or = (counts P=1) and sd = sigma
X12 = seq(0, max(h12$counts[1] + 3*sigma12, h12$counts[l12] + 3*sigma12), len=2000)
Y0_12 = dnorm(X12, h12$counts[1], sigma12)
Y1_12 = dnorm(X12, h12$counts[l12], sigma12)
# plotting dist
plot(X12, Y0_2); lines(X12, Y1_2)
# calculating BER -> area beneath the intersection of the curves * normalization
# see Fig. 1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC517707/figure/F1/
BER = sum(pmin(Y0_12, Y1_12)) * (X12[2]-X12[1])

# BER rhlE12: 8.112538e-05
