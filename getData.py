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
    
    # must be changed for future targets
    tardf = pd.read_csv("targets.csv")
    targets = [tuple(r) for r in tardf.loc[tardf['type']=='cds'][['full','short']].to_numpy()]
    extractGene.main(targets, dir="refFiles")

    # must be changed for future targets
    rnaTargets = [tuple(r) for r in tardf.loc[tardf['type']=='rna'][['full','short']].to_numpy()]
    extractGene.main(rnaTargets, dir="refFiles", type='rna')

    # must be changed for future targets
    targets = tardf['short'].to_list()
    seqClean.main(targets)


if __name__ == "__main__":
    main()
