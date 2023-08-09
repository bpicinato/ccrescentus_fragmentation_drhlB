#!/bin/bash

#####################
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
prefixes=`ls frtc_iscR/frtc_original/raw/*.fastq.gz | sed 's/frtc_iscR\/frtc_original\/raw\///;s/_R[12].*$//' | sort | uniq`
# dir with the input .bam alignment files 
coveragedir="frtc_iscR/frtc_original/coverage_uniq"
# output dir for fiveprime profile files
fiveprimedir="input_undersampling/iscR/fiveprime"
# output dir for threeprime profile files
threeprimedir="input_undersampling/iscR/threeprime"
# output temp dir
tempdir="input_undersampling/iscR/temp"
# output dir for final samples aligment files
sampledir="input_undersampling/iscR/sampled_reads"
# num of threads to use
threads=10
# size of sample
n=781455 # == ~868283*0.9 -> reads aligned in iscR_2 uniq lib * 90%

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
  
  # 2) sample n reads from each library i times
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

#### fiveprime ####
## iscR ##
# fwd
bedtools unionbedg -i $fiveprimedir/"iscR_1_S4-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"iscR_2_S5-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"iscR_3_S6-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"iscR_all-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"iscR_1_S4-5primeprofile-rev.bedgraph" \
$fiveprimedir/"iscR_2_S5-5primeprofile-rev.bedgraph" \
$fiveprimedir/"iscR_3_S6-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"iscR_all-5primeprofile-rev.bedgraph"
## NA1000 ##
# fwd
bedtools unionbedg -i $fiveprimedir/"NA_1_S1-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"NA_2_S2-5primeprofile-fwd.bedgraph" \
$fiveprimedir/"NA_3_S3-5primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"NA_all-5primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $fiveprimedir/"NA_1_S1-5primeprofile-rev.bedgraph" \
$fiveprimedir/"NA_2_S2-5primeprofile-rev.bedgraph" \
$fiveprimedir/"NA_3_S3-5primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $fiveprimedir/"NA_all-5primeprofile-rev.bedgraph"

#### 3prime ####
## iscR ##
# fwd
bedtools unionbedg -i $threeprimedir/"iscR_1_S4-3primeprofile-fwd.bedgraph" \
$threeprimedir/"iscR_2_S5-3primeprofile-fwd.bedgraph" \
$threeprimedir/"iscR_3_S6-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"iscR_all-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"iscR_1_S4-3primeprofile-rev.bedgraph" \
$threeprimedir/"iscR_2_S5-3primeprofile-rev.bedgraph" \
$threeprimedir/"iscR_3_S6-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"iscR_all-3primeprofile-rev.bedgraph"
## NA1000 ##
# fwd
bedtools unionbedg -i $threeprimedir/"NA_1_S1-3primeprofile-fwd.bedgraph" \
$threeprimedir/"NA_2_S2-3primeprofile-fwd.bedgraph" \
$threeprimedir/"NA_3_S3-3primeprofile-fwd.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"NA_all-3primeprofile-fwd.bedgraph"
# rev
bedtools unionbedg -i $threeprimedir/"NA_1_S1-3primeprofile-rev.bedgraph" \
$threeprimedir/"NA_2_S2-3primeprofile-rev.bedgraph" \
$threeprimedir/"NA_3_S3-3primeprofile-rev.bedgraph" | \
awk '{$4=$4+$5+$6}{print $1"\t"$2"\t"$3"\t"$4}' > $threeprimedir/"NA_all-3primeprofile-rev.bedgraph"

echo "Done all!"
