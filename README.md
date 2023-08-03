# Installation

## Main Dependencies/Version
Python packages: numpy/pandas, Bio.Align, matplotlib

R libraries: shiny, ggplot2, dplyr, reticulate, DT, 

Other tools: fastp, minimap2, samtools, kraken2

## Installation
To setup the pipeline after cloning from the repository, the following scripts must be run:

getData.py (as main)

prepIndex.py (as main) 

Optional if original targets are used, as the mmFiles.zip file contains the necessary files. If the targets are changed, this script is necessary.

Additionally, a kraken2 database must be built. It is recommended to build it using the genomeData folder under refFiles as a reference, since downloading default databases will be significantly more memory intensive and time consuming. On the VM, the database is labeled “ntmDB” and referenced as such in the server.R script. If the database name is different, the server.R script must be changed. 

# Usage

## Accessing/Using the Dashboard
The R Shiny dashboard can be run from bact-analysis-dev-bm-03 (currently on port 443, set in app.R file). The current method is to run kraken2 through the dashboard using the NTM kraken database, then choose a reference genome for the difference-based analysis. Results tables for both kraken2 analysis and difference-based analysis can be downloaded from the dashboard. Reads are aligned against full sequences for target genes, so changes to the PCR steps to get larger amplicons won’t require changes to the code.

## Interpreting Results
The table output has 4 columns: “diff”, “snp”, “gdiff”, and “gap”. “diff” is the total number of differences, whether mismatched bases or a gap in the reference that isn’t in the input sequence. “snp” is the number of base mismatches between the reference and input sequences. “gdiff” is the number of positions in the alignment where the reference sequence has a gap but the input doesn’t. “gap” is the total number of gaps of any length in the input sequence. 

Based on test data, near zero differences and 100% coverage is expected for correct identification, though kraken2 results should be considered for context. Other low difference/high coverage hits are expected, but will likely have more differences or lower coverage than the top hit. A cutoff of 5 differences for displaying the results is suggested, though a slider is provided to adjust it. 

The tool should differentiate species within complexes, though its use in differentiating subspecies hasn’t been tested yet. Tests also suggest that choosing the wrong reference genome won’t completely skew results and will not produce a false match with zero SNP distance and complete coverage. 

So far, kraken2 test results have either shown a clear choice for the reference genome (>50% of reads mapped to one species/subspecies) or a clear lack of a match in the NTM database (<25% of reads mapped to one species/subspecies). If results are inconclusive (e.g. most reads mapped to complex but not specific species), it may be necessary to run the distance analysis with multiple possible reference genomes and compare the results from all runs.  

# Further Steps

## Adding “de novo” Assembly
While the method for this is still unclear, de novo assembly of target sequences would be a good next step. Flye was selected as a good candidate for the actual assembly, but we are unsure how to separate reads by target before assembly without some reference-based method. Any de novo assembly approach should also be compared to this original approach to see if it can match the accuracy provided by kraken2. 

## Adding/Changing Targets
The pipeline may not be accurate for some species/subspecies identification. To fix this, it’s possible more targets are needed. I didn’t get around to making target selection possible through the dashboard, but it may be useful to do that in the future if it’s necessary to test different combinations of targets. For convenience, I’ve included a file called coreGenes.txt with every gene in the CDS files from RefSeq that’s shared across all NTM species in the data. 

To add a target to the analysis, the targets.csv file must be changed. For each target, there must be a line in the file specifying the full name of the target (as it appears in the RefSeq fasta entries), a shorter identifier for the target (preferably the gene name), and whether the target is a coding sequence (put “cds”) or RNA (put “rna”).

	
