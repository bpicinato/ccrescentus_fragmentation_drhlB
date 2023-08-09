#!/bin/bash

######################
# bpicinato Jul 2023 #
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
prefixes=`ls ../../frtc/raw/*.fastq.gz | sed 's/^..\/..\/frtc\/raw\///;s/_R1.*$//' | sort | uniq`
# dir with the input .bam alignment files 
coveragedir="../../frtc/coverage_uniq"
# output dir for fiveprime profile files
fiveprimedir="../../input_undersampling/fiveprime"
# output dir for threeprime profile files
threeprimedir="../../input_undersampling/threeprime"
# output temp dir
tempdir="../../input_undersampling/temp"
# output dir for final samples aligment files
sampledir="../../input_undersampling/sampled_reads"
# num of threads to use
threads=12
# size of sample
n=32264172 # == ~35849080*0.9 -> reads aligned in SRR13734267 * 90%

######################

# for each file, organize input and generate corresponding 5' and 3' profiles
for prefix in $prefixes ; do
  
  echo "Step 1: processing $prefix"
  
  # 1) 
  # 1.1) separating header to add back later (header should not be included in
  # the sampling process)
  samtools view -H $coveragedir/$prefix"-unpaired.bam" > $tempdir/$prefix"_header.txt"
  
  # 1.2) separating reads without header
  # obs.: num of reads in this file should be equal to the number of aligned 
  # reads in the corresponding library
  samtools view -@ $threads -b $coveragedir/$prefix"-unpaired.bam" > $tempdir/$prefix"_for-sampling.bam"
  
  
  ####
  echo "Step 2: sampling $prefix"
  
  # 2) sample n reads from each library
  samtools view -@ $threads $tempdir/$prefix"_for-sampling.bam" | \
  shuf -n $n > $tempdir/${prefix}"_"$n"reads-tmp.sam"
  # adding header to not raise error later
  cat $tempdir/$prefix"_header.txt" $tempdir/${prefix}"_"$n"reads-tmp.sam" > $tempdir/${prefix}"_"$n"reads-sampled.sam"
  
  # transfoming .sam in .bam
  samtools view -b -@ $threads \
  $tempdir/${prefix}"_"$n"reads-sampled.sam" > $sampledir/${prefix}"_"$n"reads-sampled.bam"

  
  ####
  echo "Step 3: making 5' and 3' profiles for $prefix"
  
  ## this part was adapted from steps 7-9 of frtc pipeline - https://github.com/alanlorenzetti/frtc ##
  
  # 3.1) fiveprime profile
  
  # fwd
  samtools view -@ $threads -h $sampledir/${prefix}"_"$n"reads-sampled.bam" | \
  grep "^@\|YT:Z:UU" | \
  samtools view -@ $threads -h -F 0x10 | \
  grep "^@\|XS:A:+" | \
  samtools view -@ $threads -h -b | \
  bedtools genomecov -5 -strand + -bga -ibam stdin > $fiveprimedir/$prefix"-5primeprofile-fwd.bedgraph"
  # rev
  samtools view -@ $threads -h $sampledir/${prefix}"_"$n"reads-sampled.bam" | \
  grep "^@\|YT:Z:UU" | \
  samtools view -@ $threads -h -f 0x10 | \
  grep "^@\|XS:A:-" | \
  samtools view -@ $threads -h -b | \
  bedtools genomecov -5 -strand - -bga -ibam stdin > $fiveprimedir/$prefix"-5primeprofile-rev.bedgraph"
  
  # 3.2) threeprime profile
  # fwd
  samtools view -@ $threads -h $sampledir/${prefix}"_"$n"reads-sampled.bam" | \
  grep "^@\|YT:Z:UU" | \
  samtools view -@ $threads -h -F 0x10 | \
  grep "^@\|XS:A:+" | \
  samtools view -@ $threads -h -b | \
  bedtools genomecov -3 -strand + -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-fwd.bedgraph"
  # rev
  samtools view -@ $threads -h $sampledir/${prefix}"_"$n"reads-sampled.bam" | \
  grep "^@\|YT:Z:UU" | \
  samtools view -@ $threads -h -f 0x10 | \
  grep "^@\|XS:A:-" | \
  samtools view -@ $threads -h -b | \
  bedtools genomecov -3 -strand - -bga -ibam stdin > $threeprimedir/$prefix"-3primeprofile-rev.bedgraph"
  
  echo "Done sampling part for $prefix! Removing temp files..."
  
  rm $tempdir/*

done

echo "Step 4: merging replicates"

#### 5prime ####
## WT ##
# fwd
bedtools unionbedg -i $fiveprimedir/"SRR13734262-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"SRR13734263-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"WT-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"SRR13734262-5primeprofile-rev.bedgraph" \
$fiveprimedir/"SRR13734263-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"WT-5primeprofile-rev.bedgraph"

## RhlE1 ###
# fwd
bedtools unionbedg -i $fiveprimedir/"SRR13734264-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"SRR13734265-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlE1-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"SRR13734264-5primeprofile-rev.bedgraph" \
$fiveprimedir/"SRR13734265-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlE1-5primeprofile-rev.bedgraph"

## RhlE2 ##
# fwd
bedtools unionbedg -i $fiveprimedir/"SRR13734266-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"SRR13734267-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlE2-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"SRR13734266-5primeprofile-rev.bedgraph" \
$fiveprimedir/"SRR13734267-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlE2-5primeprofile-rev.bedgraph"

## RhlE12 ###
# fwd
bedtools unionbedg -i $fiveprimedir/"SRR13734268-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"SRR13734269-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlE12-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"SRR13734268-5primeprofile-rev.bedgraph" \
$fiveprimedir/"SRR13734269-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"RhlE12-5primeprofile-rev.bedgraph"

#### 3prime ####
## WT ##
# fwd
bedtools unionbedg -i $threeprimedir/"SRR13734262-3primeprofile-fwd.bedgraph" \
$threeprimedir/"SRR13734263-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"WT-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"SRR13734262-3primeprofile-rev.bedgraph" \
$threeprimedir/"SRR13734263-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"WT-3primeprofile-rev.bedgraph"

## RhlE1 ###
# fwd
bedtools unionbedg -i $threeprimedir/"SRR13734264-3primeprofile-fwd.bedgraph" \
$threeprimedir/"SRR13734265-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlE1-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"SRR13734264-3primeprofile-rev.bedgraph" \
$threeprimedir/"SRR13734265-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlE1-3primeprofile-rev.bedgraph"

## RhlE2 ##
# fwd
bedtools unionbedg -i $threeprimedir/"SRR13734266-3primeprofile-fwd.bedgraph" \
$threeprimedir/"SRR13734267-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlE2-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"SRR13734266-3primeprofile-rev.bedgraph" \
$threeprimedir/"SRR13734267-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlE2-3primeprofile-rev.bedgraph"

## RhlE12 ###
# fwd
bedtools unionbedg -i $threeprimedir/"SRR13734268-3primeprofile-fwd.bedgraph" \
$threeprimedir/"SRR13734269-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlE12-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"SRR13734268-3primeprofile-rev.bedgraph" \
$threeprimedir/"SRR13734269-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"RhlE12-3primeprofile-rev.bedgraph"

echo "Done all!"
