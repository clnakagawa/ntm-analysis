import pandas as pd
import numpy as np
import cleanFTP, getGenes, getRNA, getGenomes, extractRNA, extractGene, extractITS, seqClean, os

def main():
    #clean ftp table
    print("finding relevant entries...")
    cleanFTP.main()

    #get cds, genomes, rna
    print("retrieving cds files...")
    getGenes.main()
    print("retrieving rna files...")
    getRNA.main()
    print("retrieving genome files...")
    getGenomes.main()

    # define targets
    gene_targets = [('chaperonin GroEL', 'groL'),
               ('ATP synthase subunit beta', 'atpD'),
               ('RNA polymerase subunit beta', 'rpoB'),
               ('elongation factor Tu', 'tuf')]
    rna_targets = [('16S', '16S'),('23S','23S')]

    # extract targets from files
    print("finding targets in files...")
    extractGene.main(gene_targets)
    extractRNA.main(rna_targets)
    print("finding ITS...")
    extractITS.main()

    # clean files (remove dupes/bad matches)
    print("cleaning target lists...")
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf']
    seqClean.main(targets)
    print("data ready for alignment")
    print("launching alignment")
    os.chdir("C:\\Users\cyn06\OneDrive - New York State Office of Information Technology Services\Downloads\mafft-7.520-win64-signed\mafft-win")
    os.startfile('allAlign.bat')

if __name__ == "__main__":
    main()
