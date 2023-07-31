#!/bin/bash
mkdir fqProcessed/$4

fastp -i $1 -A > -o $2/trimmed.fastq

minimap2 -a "mmInds/$2/$3_16S.mmi" $1 > fqProcessed/$4/16Salign.sam
samtools sort fqProcessed/$4/16Salign.sam -o fqProcessed/$4/16Salign.bam
samtools consensus -f fasta fqProcessed/$4/16Salign.bam -o fqProcessed/$4/16Sconsensus.fasta
samtools depth fqProcessed/$4/16Salign.bam > fqProcessed/$4/16Scoverage.txt

minimap2 -a "mmInds/$2/$3_23S.mmi" $1 > fqProcessed/$4/23Salign.sam
samtools sort fqProcessed/$4/23Salign.sam -o fqProcessed/$4/23Salign.bam
samtools consensus -f fasta fqProcessed/$4/23Salign.bam -o fqProcessed/$4/23Sconsensus.fasta
samtools depth fqProcessed/$4/23Salign.bam > fqProcessed/$4/23Scoverage.txt

minimap2 -a "mmInds/$2/$3_atpD.mmi" $1 > fqProcessed/$4/atpDalign.sam
samtools sort fqProcessed/$4/atpDalign.sam -o fqProcessed/$4/atpDalign.bam
samtools consensus -f fasta fqProcessed/$4/atpDalign.bam -o fqProcessed/$4/atpDconsensus.fasta
samtools depth fqProcessed/$4/atpDalign.bam > fqProcessed/$4/atpDcoverage.txt

minimap2 -a "mmInds/$2/$3_groL.mmi" $1 > fqProcessed/$4/groLalign.sam
samtools sort fqProcessed/$4/groLalign.sam -o fqProcessed/$4/groLalign.bam
samtools consensus -f fasta fqProcessed/$4/groLalign.bam -o fqProcessed/$4/groLconsensus.fasta
samtools depth fqProcessed/$4/groLalign.bam > fqProcessed/$4/groLcoverage.txt

minimap2 -a "mmInds/$2/$3_rpoB.mmi" $1 > fqProcessed/$4/rpoBalign.sam
samtools sort fqProcessed/$4/rpoBalign.sam -o fqProcessed/$4/rpoBalign.bam
samtools consensus -f fasta fqProcessed/$4/rpoBalign.bam -o fqProcessed/$4/rpoBconsensus.fasta
samtools depth fqProcessed/$4/rpoBalign.bam > fqProcessed/$4/rpoBcoverage.txt

minimap2 -a "mmInds/$2/$3_tuf.mmi" $1 > fqProcessed/$4/tufalign.sam
samtools sort fqProcessed/$4/tufalign.sam -o fqProcessed/$4/tufalign.bam
samtools consensus -f fasta fqProcessed/$4/tufalign.bam -o fqProcessed/$4/tufconsensus.fasta
samtools depth fqProcessed/$4/tufalign.bam > fqProcessed/$4/tufcoverage.txt

cat fqProcessed/$4/16Sconsensus.fasta fqProcessed/$4/23Sconsensus.fasta fqProcessed/$4/atpDconsensus.fasta fqProcessed/$4/groLconsensus.fasta fqProcessed/$4/rpoBconsensus.fasta fqProcessed/$4/tufconsensus.fasta > fqProcessed/$4/allSeqs.txt

