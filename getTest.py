import pandas as pd
import numpy as np
import cleanFTP, getGenes, getRNA, getGenomes, random, seqClean
import extractGene
import extractRNA


def main():
    # cleanFTP.test()
    # getGenes.test()
    # getRNA.test()
    # getGenomes.test()

    # check if all targets present
    # get accession numbers
    # OLD CODE from getting the list of ones w targets
    # df = pd.read_csv("testFiles/summaryData.csv")
    # accnums = df['# assembly_accession']
    # geneNames = ['ATP synthase subunit beta', 'chaperonin GroEL', 'RNA polymerase subunit beta', 'elongation factor Tu']
    # rnaNames = ['16S', '23S']
    #
    # for target in geneNames:
    #     print(target)
    #     acctemp = []
    #     for a in accnums:
    #         with open(f"testFiles/cdsData/{a}_cds.fna", 'r') as f:
    #             if target in f.read():
    #                 acctemp.append(a)
    #     accnums = acctemp
    # for target in rnaNames:
    #     print(target)
    #     acctemp = []
    #     for a in accnums:
    #         with open(f"testFiles/rnaData/{a}_rna.fna", 'r') as f:
    #             if target in f.read():
    #                 acctemp.append(a)
    #     accnums = acctemp
    #
    # print(len(accnums))
    # df = pd.DataFrame({'accs':accnums})
    # df.to_csv("testFiles/hasTargetList.csv")

    df = pd.read_csv('testFiles/hasTargetList.csv')
    accList = list(df['accs'])

    # get random sample
    toGet = random.sample(accList, 20)
    print(toGet)
    with open("testFiles/testList.txt", 'w') as f:
        f.write('\n'.join(toGet))

    # pull sequences and send to query
    genetargets = [('chaperonin GroEL', 'groL'),
               ('ATP synthase subunit beta', 'atpD'),
               ('RNA polymerase subunit beta', 'rpoB'),
               ('elongation factor Tu', 'tuf')]
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf']

    extractGene.test(genetargets, toGet)
    rnatargets = [('16S', '16S'), ('23S', '23S')]
    extractRNA.test(rnatargets, toGet)
    seqClean.main(targets, test=True)

if __name__ == "__main__":
    main()