# ntm-analysis

Installation

Main Dependencies/Version
Python packages: numpy/pandas, Bio.Align, matplotlib
R libraries: shiny, ggplot2, dplyr, reticulate, DT, 
Other tools: fastp, minimap2, samtools, kraken2

Installation
	The pipeline relies on data pulled from RefSeq, used to build databases/datasets and as a reference itself. The getData.py file pulls all the necessary files for the python scripts. Additionally, a kraken2 database must be built, for which the genomeData folder containing all NTM genomes downloaded from RefSeq should be used. For building minimap indexes, I have a folder called mmRef with a folder for every species containing a file for each target sequence.

Usage

Accessing/Using the Dashboard
The R Shiny dashboard can be run from bact-analysis-dev-bm-03 (currently on port 443, set in app.R file). The current procedure is to run kraken2 through the dashboard using the NTM krakenDB, then choose a reference genome for the difference-based analysis. Results tables for both kraken2 analysis and difference-based analysis can be downloaded from the dashboard. Reads are aligned against full sequences for target genes, so changes to the PCR steps to get larger amplicons won’t require changes to the code.

Interpreting Results

Based on test data, near zero SNP distance and 100% coverage is expected for correct identification, though kraken2 results should be considered for context. Other low SNP/high coverage hits are expected, but will likely have higher SNP distance or lower coverage than the top hit. A cutoff of 5 SNP distance for displaying the results is suggested, though a slider is provided to adjust it. The number of gap positions for the amplified region is also listed and should be considered along with the SNP distance for interpreting results. 

The tool should differentiate species within complexes, though its use in differentiating subspecies hasn’t been tested yet. Tests also suggest that choosing the wrong reference genome won’t completely skew results and will not produce a false match with zero SNP distance and complete coverage. Kraken2 test results have either shown a clear choice for the reference genome (>50% of reads mapped to one species/subspecies) or a clear lack of a match in the NTM database (<25% of reads mapped to one species/subspecies). If results are inconclusive (e.g. most reads mapped to complex but not specific species), it may be necessary to run the distance analysis with multiple possible reference genomes and compare the results from all runs.  

Further Steps

Adding “de novo” Assembly

The method for this is still unclear, but de novo assembly of target sequences would be a good next step. Flye was selected as a good candidate for the actual assembly, but we are unsure how to separate reads by target before assembly without some reference-based method. Any de novo assembly approach should also be compared to this original approach to see if it can match the accuracy provided by kraken2. 

Adding/Changing Targets

The pipeline may not be accurate for many species/subspecies identification. To fix this, it’s possible more targets are needed. While I unfortunately didn’t get around to making target selection possible through the dashboard, it may be useful to do that in the future if it’s necessary to test different combinations of targets. To get a list of possible targets, I’ve included a file called coreGenes.txt with every gene in the CDS files from RefSeq that’s shared across all NTM species in the data. 

To add a target to the analysis, the targets.csv file must be changed. For each target, there must be a line in the file specifying the full name of the target (as it appears in the RefSeq fasta entries), a shorter identifier for the target (preferably the gene name), and whether the target is a coding sequence (put “cds”) or RNA (put “rna”).

	
