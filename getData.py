import pandas as pd
import numpy as np
import cleanFTP, getGenes, getRNA, getGenomes, extractRNA, extractGene, extractITS, seqClean, os

def main(file=None):
    #clean ftp table
    print("finding relevant entries...")
    if file is None:
        cleanFTP.main()

    #get cds, genomes, rna
    print("retrieving files...")
    getGenes.main(type="cds")
    getGenes.main(type="rna")
    getGenes.main()

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

if __name__ == "__main__":
    main()
