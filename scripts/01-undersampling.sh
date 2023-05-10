#!/bin/bash

######################
# bpicinato Apr 2023 #
######################

######################
#   DESCRIPTION     #
######################
# script to generate a sample n reads from a .bam alignment file and then make
# the correspondent 5' and 3' profiles as in frtc alingment pipeline
# https://github.com/alanlorenzetti/frtc

# used to make undersampling from different RNA-seq libraries so they could be
# compared in degradation analysis

######################

# getting file names and prefixes
prefixes=`ls frtc_drhlB/ccrescentus_hugo_v2/raw/*.fastq.gz | sed 's/^frtc_drhlB\/ccrescentus_hugo_v2\/raw\///;s/_R[12].*$//' | sort | uniq`
# dir with the input .bam alignment files 
coveragedir="frtc_drhlB/ccrescentus_hugo_v2/coverage_uniq"
# output dir for fiveprime profile files
fiveprimedir="degradation_rhlB/input_undersampling/fiveprime"
# output dir for threeprime profile files
threeprimedir="degradation_rhlB/input_undersampling/threeprime"
# output temp dir
tempdir="degradation_rhlB/input_undersampling/temp"
# output dir for final samples aligment files
sampledir="degradation_rhlB/input_undersampling/sampled_reads"
# num of threads to use
threads=10
# size of sample
n=997243 # == ~1108048*0.9 -> reads aligned in RhlB1_30C lib * 90%

######################

# for each file, organize input and generate corresponding 5' and 3' profiles
for prefix in $prefixes ; do
  
  echo "Step 1: processing $prefix"
  
  # 1) merge paired and unpaired BAM files from coverage folder to make one big
  #    file containing all mapped reads without multimappers (uniq folder)
  samtools merge -@ $threads \
  $tempdir/$prefix"-all.bam" \
  $coveragedir/$prefix"-paired.bam" \
  $coveragedir/$prefix"-unpaired.bam"

  # 1.1) separating header to add back later (header should not be included in
  # the sampling process)
  samtools view -H $tempdir/$prefix"-all.bam" > $tempdir/$prefix"_header.txt"
  
  # 1.2) separating reads without R2 paired - for sampling only one per pair
  # obs.: num of reads in this file should be equal to the number of aligned 
  # reads in the corresponding library
  samtools view -@ $threads -b -F 0x80 $tempdir/$prefix"-all.bam" > $tempdir/$prefix"_for-sampling.bam"
  
  # 1.3) separating paired R2 reads to be added later, after sampling
  samtools view -@ $threads -b -f 0x80 $tempdir/$prefix"-all.bam" > $tempdir/$prefix"-paired-R2.bam"
  
  ####
  echo "Step 2: sampling $prefix"
  
  # 2) sample n reads from each library
  samtools view -@ $threads $tempdir/$prefix"_for-sampling.bam" | \
  shuf -n $n > $tempdir/${prefix}"_"$n"reads-R1-tmp.sam"
  # adding header to not raise error later
  cat $tempdir/$prefix"_header.txt" $tempdir/${prefix}"_"$n"reads-R1-tmp.sam" > $tempdir/${prefix}"_"$n"reads-R1.sam"
  
  # getting R1 read names to get pair in R2 file
  samtools view -@ $threads -f 0x40 $tempdir/${prefix}"_"$n"reads-R1.sam" | \
  awk {'print $1'} > $tempdir/$prefix"_reads-to-get.txt"
  # grepping corresponding lines in R2 only file
  samtools view -@ $threads $tempdir/$prefix"-paired-R2.bam" | \
  grep -f $tempdir/$prefix"_reads-to-get.txt" > $tempdir/${prefix}"_"$n"reads-R2-tmp.sam"
  # adding header to not raise error later
  cat $tempdir/$prefix"_header.txt" $tempdir/${prefix}"_"$n"reads-R2-tmp.sam" > $tempdir/${prefix}"_"$n"reads-R2.sam"
  
  # merging R1 and R2 in one file
  samtools merge -@ $threads \
  $sampledir/${prefix}"_"$n"reads_all.bam" \
  $tempdir/${prefix}"_"$n"reads-R1.sam" \
  $tempdir/${prefix}"_"$n"reads-R2.sam"
  
  ####
  echo "Step 3: making 5' and 3' profiles for $prefix"
  
  ## this part was adapted from steps 7-9 of frtc pipeline - https://github.com/alanlorenzetti/frtc ##
  
  # 3.1) fiveprime profile
  # generating profile for paired R1 (flag 0x40 -> first in pair, only in paired reads)
  # fwd
  samtools view -@ $threads -h -b -f 0x40 $sampledir/${prefix}"_"$n"reads_all.bam" | \
  bedtools genomecov -5 -strand + -bga -ibam stdin > $fiveprimedir/$prefix"-5primeprofile-paired-fwd.bedgraph"
  # rev
  samtools view -@ $threads -h -b -f 0x40 $sampledir/${prefix}"_"$n"reads_all.bam" | \
  bedtools genomecov -5 -strand - -bga -ibam stdin > $fiveprimedir/$prefix"-5primeprofile-paired-rev.bedgraph"
  
  # generating profile for unpaired reads
  # fwd
  samtools view -@ $threads -h $sampledir/${prefix}"_"$n"reads_all.bam" | \
  grep "^@\|YT:Z:UU" | \
  samtools view -@ $threads -h -F 0x10 | \
  awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 !~ /_R2$/){print}}}' | \
  grep "^@\|XS:A:+" | \
  samtools view -@ $threads -h -b | \
  bedtools genomecov -5 -strand + -bga -ibam stdin > $fiveprimedir/$prefix"-5primeprofile-unpaired-fwd.bedgraph"
  # rev
  samtools view -@ $threads -h $sampledir/${prefix}"_"$n"reads_all.bam" | \
  grep "^@\|YT:Z:UU" | \
  samtools view -@ $threads -h -f 0x10 | \
  awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 !~ /_R2$/){print}}}' | \
  grep "^@\|XS:A:-" | \
  samtools view -@ $threads -h -b | \
  bedtools genomecov -5 -strand - -bga -ibam stdin > $fiveprimedir/$prefix"-5primeprofile-unpaired-rev.bedgraph"
  
  # merging paired and unpaired
  bedtools unionbedg -i \
  $fiveprimedir/$prefix"-5primeprofile-paired-fwd.bedgraph" \
  $fiveprimedir/$prefix"-5primeprofile-unpaired-fwd.bedgraph" | \
  awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,($4+$5)}' > $fiveprimedir/$prefix"-5primeprofile-fwd.bedgraph"

  bedtools unionbedg -i \
  $fiveprimedir/$prefix"-5primeprofile-paired-rev.bedgraph" \
  $fiveprimedir/$prefix"-5primeprofile-unpaired-rev.bedgraph" | \
  awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,($4+$5)}' > $fiveprimedir/$prefix"-5primeprofile-rev.bedgraph"
  
  # 3.2) threeprime profile
  # generating profile for paired R2 reads (flag 0x80 -> second in pair, only in paired reads)
  # fwd
  samtools view -@ $threads -h -b -f 0x80 $sampledir/${prefix}"_"$n"reads_all.bam" | \
  bedtools genomecov -5 -strand - -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-paired-fwd.bedgraph"
  # rev
  samtools view -@ $threads -h -b -f 0x80 $sampledir/${prefix}"_"$n"reads_all.bam" | \
  bedtools genomecov -5 -strand + -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-paired-rev.bedgraph"
  
  # generating profile for unpaired reads
  # fwd
  samtools view -@ $threads -h $sampledir/${prefix}"_"$n"reads_all.bam" | \
  grep "^@\|YT:Z:UU" | \
  samtools view -@ $threads -h -F 0x10 | \
  awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 ~ /_R2$/){print}}}' | grep "^@\|XS:A:+" | \
  samtools view -@ $threads -h -b | \
  bedtools genomecov -3 -strand + -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-unpaired-fwd.bedgraph"
  # rev
  samtools view -@ $threads -h $sampledir/${prefix}"_"$n"reads_all.bam" | \
  grep "^@\|YT:Z:UU" | \
  samtools view -@ $threads -h -f 0x10 | \
  awk -v OFS="\t" -v FS="\t" '{if(/^@/){print}else{if($1 ~ /_R2$/){print}}}' | grep "^@\|XS:A:-" | \
  samtools view -@ $threads -h -b | \
  bedtools genomecov -3 -strand - -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-unpaired-rev.bedgraph"

  # merging paired and unpaired
  bedtools unionbedg -i \
  $threeprimedir/$prefix"-3primeprofile-paired-fwd.bedgraph" \
  $threeprimedir/$prefix"-3primeprofile-unpaired-fwd.bedgraph" | \
  awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,($4+$5)}' > $threeprimedir/$prefix"-3primeprofile-fwd.bedgraph"

  bedtools unionbedg -i \
  $threeprimedir/$prefix"-3primeprofile-paired-rev.bedgraph" \
  $threeprimedir/$prefix"-3primeprofile-unpaired-rev.bedgraph" | \
  awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,($4+$5)}' > $threeprimedir/$prefix"-3primeprofile-rev.bedgraph"
  
  echo "Done sampling part for $prefix! Removing temp files..."
  
  rm $tempdir/*

done

echo "Step 4: merging replicates"

#### 5prime ####
### 10C ###
## RhlB ##
# fwd
bedtools unionbedg -i $fiveprimedir/"RhlB1_10C-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"RhlB2_10C-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"RhlB3_10C-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlB_all_10C-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"RhlB1_10C-5primeprofile-rev.bedgraph" \
$fiveprimedir/"RhlB2_10C-5primeprofile-rev.bedgraph" \
$fiveprimedir/"RhlB3_10C-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlB_all_10C-5primeprofile-rev.bedgraph"
## NA1000 ###
# fwd
bedtools unionbedg -i $fiveprimedir/"NA1_10C-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"NA2_10C-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"NA3_10C-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"NA_all_10C-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"NA1_10C-5primeprofile-rev.bedgraph" \
$fiveprimedir/"NA2_10C-5primeprofile-rev.bedgraph" \
$fiveprimedir/"NA3_10C-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"NA_all_10C-5primeprofile-rev.bedgraph"

### 30C ###
## RhlB ##
# fwd
bedtools unionbedg -i $fiveprimedir/"RhlB1_30C-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"RhlB2_30C-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"RhlB3_30C-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlB_all_30C-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"RhlB1_30C-5primeprofile-rev.bedgraph" \
$fiveprimedir/"RhlB2_30C-5primeprofile-rev.bedgraph" \
$fiveprimedir/"RhlB3_30C-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlB_all_30C-5primeprofile-rev.bedgraph"
## NA1000 ###
# fwd
bedtools unionbedg -i $fiveprimedir/"NA1_30C-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"NA2_30C-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"NA3_30C-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"NA_all_30C-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"NA1_30C-5primeprofile-rev.bedgraph" \
$fiveprimedir/"NA2_30C-5primeprofile-rev.bedgraph" \
$fiveprimedir/"NA3_30C-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"NA_all_30C-5primeprofile-rev.bedgraph"

#### 3prime ####
### 10C ###
## RhlB ##
# fwd
bedtools unionbedg -i $threeprimedir/"RhlB1_10C-3primeprofile-fwd.bedgraph" \
$threeprimedir/"RhlB2_10C-3primeprofile-fwd.bedgraph" \
$threeprimedir/"RhlB3_10C-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlB_all_10C-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"RhlB1_10C-3primeprofile-rev.bedgraph" \
$threeprimedir/"RhlB2_10C-3primeprofile-rev.bedgraph" \
$threeprimedir/"RhlB3_10C-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlB_all_10C-3primeprofile-rev.bedgraph"
## NA1000 ##
# fwd
bedtools unionbedg -i $threeprimedir/"NA1_10C-3primeprofile-fwd.bedgraph" \
$threeprimedir/"NA2_10C-3primeprofile-fwd.bedgraph" \
$threeprimedir/"NA3_10C-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"NA_all_10C-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"NA1_10C-3primeprofile-rev.bedgraph" \
$threeprimedir/"NA2_10C-3primeprofile-rev.bedgraph" \
$threeprimedir/"NA3_10C-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"NA_all_10C-3primeprofile-rev.bedgraph"

### 30C ###
## RhlB ##
# fwd
bedtools unionbedg -i $threeprimedir/"RhlB1_30C-3primeprofile-fwd.bedgraph" \
$threeprimedir/"RhlB2_30C-3primeprofile-fwd.bedgraph" \
$threeprimedir/"RhlB3_30C-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlB_all_30C-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"RhlB1_30C-3primeprofile-rev.bedgraph" \
$threeprimedir/"RhlB2_30C-3primeprofile-rev.bedgraph" \
$threeprimedir/"RhlB3_30C-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlB_all_30C-3primeprofile-rev.bedgraph"
## NA1000 ##
# fwd
bedtools unionbedg -i $threeprimedir/"NA1_30C-3primeprofile-fwd.bedgraph" \
$threeprimedir/"NA2_30C-3primeprofile-fwd.bedgraph" \
$threeprimedir/"NA3_30C-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"NA_all_30C-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"NA1_30C-3primeprofile-rev.bedgraph" \
$threeprimedir/"NA2_30C-3primeprofile-rev.bedgraph" \
$threeprimedir/"NA3_30C-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"NA_all_30C-3primeprofile-rev.bedgraph"

echo "Done all!"
