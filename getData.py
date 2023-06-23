import pandas as pd
import numpy as np
import cleanFTP, getGenes, extractGene, extractITS, seqClean, os

def main(file=None):
    #clean ftp table
    print("finding relevant entries...")
    if file is None:
        cleanFTP.main()

    #get cds, genomes, rna
    print("retrieving files...")
    getGenes.main(output="refFiles",type="cds")
    getGenes.main(output="refFiles",type="rna")
    getGenes.main(output="refFiles")

    targets = [('chaperonin GroEL', 'groL'),
               ('ATP synthase subunit beta', 'atpD'),
               ('RNA polymerase subunit beta', 'rpoB'),
               ('elongation factor Tu', 'tuf')]
    extractGene.main(targets, dir="refFiles")
    rnaTargets = [('16S','16S'),('23S','23S')]
    extractGene.main(rnaTargets, dir="refFiles", type='rna')

    targets = ['groL', 'atpD', 'rpoB', 'tuf', '16S', '23S']
    seqClean.main(targets)


if __name__ == "__main__":
    main()
