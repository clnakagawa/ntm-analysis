#!/bin/bash
mkdir $4

export BOWTIE2_INDEXES=~/Users/cln/BTInds
fastp -i $1 > -o $2/trimmed.fastq

./minimap2/minimap2 -a "mmInds/$2/$3_16S.mmi" $1 > $4/16Salign.sam
samtools sort $4/16Salign.sam -o $4/16Salign.bam
samtools consensus -f fasta $4/16Salign.bam -o $4/16Sconsensus.fasta

./minimap2/minimap2 -a "mmInds/$2/$3_23S.mmi" $1 > $4/23Salign.sam
samtools sort $4/23Salign.sam -o $4/23Salign.bam
samtools consensus -f fasta $4/23Salign.bam -o $4/23Sconsensus.fasta

./minimap2/minimap2 -a "mmInds/$2/$3_atpD.mmi" $1 > $4/atpDalign.sam
samtools sort $4/atpDalign.sam -o $4/atpDalign.bam
samtools consensus -f fasta $4/atpDalign.bam -o $4/atpDconsensus.fasta

./minimap2/minimap2 -a "mmInds/$2/$3_groL.mmi" $1 > $4/groLalign.sam
samtools sort $4/groLalign.sam -o $4/groLalign.bam
samtools consensus -f fasta $4/groLalign.bam -o $4/groLconsensus.fasta

./minimap2/minimap2 -a "mmInds/$2/$3_rpoB.mmi" $1 > $4/rpoBalign.sam
samtools sort $4/rpoBalign.sam -o $4/rpoBalign.bam
samtools consensus -f fasta $4/rpoBalign.bam -o $4/rpoBconsensus.fasta

./minimap2/minimap2 -a "mmInds/$2/$3_tuf.mmi" $1 > $4/tufalign.sam
samtools sort $4/tufalign.sam -o $4/tufalign.bam
samtools consensus -f fasta $4/tufalign.bam -o $4/tufconsensus.fasta

cat $4/16Sconsensus.fasta $4/23Sconsensus.fasta $4/atpDconsensus.fasta $4/groLconsensus.fasta $4/rpoBconsensus.fasta $4/tufconsensus.fasta > $4/allSeqs.txt

