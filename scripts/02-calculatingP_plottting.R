######################
# bpicinato Apr 2023 #
######################

######################
#   DESCRIPTION     #
######################
# script to calculate values of P for whole genome and generate plot

######################

library("tidyverse")
library("ggplot2")

options(scipen=10000)

# sample 1 used for plot

# ATTENTION: using "threeprime" files that are actually fiveprime of transcript
# because of sequencing kit and protocol - dUTP protocol
# (Illumina Stranded Total RNA Prep Ligation with Ribozero Plus)
# also, strandness is inverted

###########
### 10C ###
###########

# reading

##### inverted because of kit - 3 prime files, but actually real 5 prime #####

NA10C_fwd_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/rhlB/uniq/sample1/threeprime/NA_all_10C-3primeprofile-fwd.bedgraph", col_names = F)
NA10C_rev_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/rhlB/uniq/sample1/threeprime/NA_all_10C-3primeprofile-rev.bedgraph", col_names = F)
RhlB10C_fwd_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/rhlB/uniq/sample1/threeprime/RhlB_all_10C-3primeprofile-fwd.bedgraph", col_names = F)
RhlB10C_rev_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/rhlB/uniq/sample1/threeprime/RhlB_all_10C-3primeprofile-rev.bedgraph", col_names = F)

NA10C_fwd_3 = rep(c(NA10C_fwd_bedg_3$X4), times = NA10C_fwd_bedg_3$X3 - NA10C_fwd_bedg_3$X2)
NA10C_rev_3 = rep(c(NA10C_rev_bedg_3$X4), times = NA10C_rev_bedg_3$X3 - NA10C_rev_bedg_3$X2)

RhlB10C_fwd_3 = rep(c(RhlB10C_fwd_bedg_3$X4), times = RhlB10C_fwd_bedg_3$X3 - RhlB10C_fwd_bedg_3$X2)
RhlB10C_rev_3 = rep(c(RhlB10C_rev_bedg_3$X4), times = RhlB10C_rev_bedg_3$X3 - RhlB10C_rev_bedg_3$X2)

p_fwd_3_10C = RhlB10C_fwd_3 / (RhlB10C_fwd_3 + NA10C_fwd_3)
p_rev_3_10C = RhlB10C_rev_3 / (RhlB10C_rev_3 + NA10C_rev_3)
p_threeprime_10C = c(p_fwd_3_10C, p_rev_3_10C)

###########
### 30C ###
###########

# reading

NA30C_fwd_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/rhlB/uniq/sample1/threeprime/NA_all_30C-3primeprofile-fwd.bedgraph", col_names = F)
NA30C_rev_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/rhlB/uniq/sample1/threeprime/NA_all_30C-3primeprofile-rev.bedgraph", col_names = F)
RhlB30C_fwd_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/rhlB/uniq/sample1/threeprime/RhlB_all_30C-3primeprofile-fwd.bedgraph", col_names = F)
RhlB30C_rev_bedg_3 = read_tsv("degradation_rhlB/input_undersampling/rhlB/uniq/sample1/threeprime/RhlB_all_30C-3primeprofile-rev.bedgraph", col_names = F)

NA30C_fwd_3 = rep(c(NA30C_fwd_bedg_3$X4), times = NA30C_fwd_bedg_3$X3 - NA30C_fwd_bedg_3$X2)
NA30C_rev_3 = rep(c(NA30C_rev_bedg_3$X4), times = NA30C_rev_bedg_3$X3 - NA30C_rev_bedg_3$X2)

RhlB30C_fwd_3 = rep(c(RhlB30C_fwd_bedg_3$X4), times = RhlB30C_fwd_bedg_3$X3 - RhlB30C_fwd_bedg_3$X2)
RhlB30C_rev_3 = rep(c(RhlB30C_rev_bedg_3$X4), times = RhlB30C_rev_bedg_3$X3 - RhlB30C_rev_bedg_3$X2)

p_fwd_3_30C = RhlB30C_fwd_3 / (RhlB30C_fwd_3 + NA30C_fwd_3)
p_rev_3_30C = RhlB30C_rev_3 / (RhlB30C_rev_3 + NA30C_rev_3)
p_threeprime_30C = c(p_fwd_3_30C, p_rev_3_30C)

#######

# putting into dataframe to plot
df = data.frame("10C" = p_threeprime_10C, "30C" = p_threeprime_30C)
df_pivot = df %>% 
  pivot_longer(everything(), names_to = "temperature", values_to = "P" )

# making dataframe with horizontal line information for plot
lines = data.frame(temperature = c("X10C", "X30C"), 
                   line = c(sum(abs(p_threeprime_10C - 0) < 1e-6, na.rm = T), sum(abs(p_threeprime_30C - 0) < 1e-6, na.rm = T)))
# labels for facets 
# https://community.rstudio.com/t/subscripts-and-superscripts-facet-wrap-facet-labels/81072/6
labels = as_labeller(c(X10C = "10~degree*C", 
                       X30C = "30~degree*C"), default = label_parsed)
# plotting
ggplot(df_pivot, aes(x = P, fill = factor(temperature))) + 
  geom_histogram(bins = 20, color = "white") + 
  facet_wrap(~temperature, labeller = labels) +
  geom_abline(data = lines, aes(intercept = line, slope = 0), linetype = "dashed", 
              color = "grey", size = 0.25) + 
  theme_minimal() +
  scale_fill_manual(values = c("#4471c4", "#eb7c3b")) +
  labs(x="P", y = "Frequency") +
  theme(panel.border = element_rect(linetype = "solid", color = "black", fill = NA),
        legend.position = "none",
        text = element_text(color = "black", size = 8),
        axis.text = element_text(color = "black"), 
        strip.text.x = element_text(size = 10))
# saving
ggsave("degradation_rhlB/final_plots/P-5prime-profile.png", 
       dpi = 600, width = 6, height = 3, units = "in", bg = "white")
  

 
