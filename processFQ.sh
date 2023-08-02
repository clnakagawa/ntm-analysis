#!/bin/bash
mkdir fqProcessed/$4

fastp -i $1 -A > -o $2/trimmed.fastq

FNAME=$1
SP_NAME=$2
IND=$3
SAMPLE_DIR=$4

fastp -i $FNAME -A > -o $SP_NAME/trimmed.fastq

echo $IND_NAME

map_to_consensus()
{
   minimap2 -a mmInds/$SP_NAME/${IND}_$1.mmi $FNAME > fqProcessed/$SAMPLE_DIR/$1align.sam
   samtools sort fqProcessed/$SAMPLE_DIR/$1align.sam -o fqProcessed/$SAMPLE_DIR/$1align.bam
   samtools consensus -f fasta fqProcessed/$SAMPLE_DIR/$1align.bam -o fqProcessed/$SAMPLE_DIR/$1consensus.fasta
   samtools depth fqProcessed/$SAMPLE_DIR/$1align.bam > fqProcessed/$SAMPLE_DIR/$1coverage.txt
}

while IFS="," read -r rec1
do 
   map_to_consensus $rec1
done < <(cut -d "," -f2 targets.csv | tail -n +2)

cat fqProcessed/$SAMPLE_DIR/*consensus.fasta > fqProcessed/$SAMPLE_DIR/allSeqs.txt
