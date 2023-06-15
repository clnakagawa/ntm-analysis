import pandas as pd
from Bio import Align
import cleanFTP, getGenes, getRNA, getGenomes, random, seqClean, os, snpDists, datetime
import extractGene
import extractRNA
from functools import reduce

def main():
    pd.set_option('display.max_columns', 15)
    # cleanFTP.test()
    # getGenes.test()
    # getRNA.test()
    # getGenomes.test()

    # check if all targets present
    # get accession numbers
    # OLD CODE from getting the list of ones w targets
    sumData = pd.read_csv("testFiles/summaryData.csv")
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
    try:
        with open("query/toTest.txt", 'r') as f:
            accList = f.read().split("\n")
    except:
        accList = list(df['accs'])

    try:
        resdf = pd.read_csv('query/queryResults.csv', index_col=0)
    except:
        resdf = pd.DataFrame(columns=['acc', 'name', 'id', 'realDist', 'idDist'])


    while accList:
        acc = accList.pop()
        start = datetime.datetime.now()
        if not os.path.exists(f"query/{acc}"):
            os.mkdir(f"query/{acc}")
        # pull sequences and send to query
        genetargets = [('chaperonin GroEL', 'groL'),
                   ('ATP synthase subunit beta', 'atpD'),
                   ('RNA polymerase subunit beta', 'rpoB'),
                   ('elongation factor Tu', 'tuf')]
        targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf']

        extractGene.test(genetargets, acc)
        rnatargets = [('16S', '16S'), ('23S', '23S')]
        extractRNA.test(rnatargets, acc)
        seqClean.test(targets, acc)
        dfs = []
        aligner = Align.PairwiseAligner()
        if not os.path.exists(f"queryAligns/{acc}"):
            os.mkdir(f"queryAligns/{acc}")
        for target in targets:
            with open(f"query/{acc}/{target}.txt") as f:
                qseq = ''.join(f.read().split("\n")[1:])
            with open(f"cleanedSeqs/new{target}_allseqs.txt", 'r') as f:
                refs = f.read().split("\n>")
            seqs = [''.join(seq.split("\n")[1:]) for seq in refs]
            sps = [seq.split("\n")[0].replace(">","") for seq in refs]
            alns = [aligner.align(qseq, seq)[0] for seq in seqs]
            snps = snpDists.distanceList([aln[0] for aln in alns],
                                         [aln[1] for aln in alns])
            df = pd.DataFrame({'sp':sps,'snp':snps}).sort_values('snp')
            df = df.set_index(['sp'])
            df.to_csv(f"query/{acc}/{target}SNP.csv")
            dfs.append(df)
        total = reduce(lambda x, y: x.add(y), dfs)
        total = total.sort_values('snp')
        total.to_csv(f"query/{acc}/totalSNP.csv")
        print(total.head())
        realEntry = sumData.loc[sumData['# assembly_accession'] == acc]
        print(acc,
              realEntry['organism_name'].values[0],
              total.index.values[0],
              total.loc[realEntry['sp'].values[0]]['snp'],
              total['snp'].values[0])
        resdf.loc[len(resdf.index)] = [acc,
                                       realEntry['organism_name'].values[0],
                                       total.index.values[0],
                                       total.loc[realEntry['sp'].values[0]]['snp'],
                                       total['snp'].values[0]]

        print(f"Time is {datetime.datetime.now()-start}")
        with open("query/toTest.txt", 'w') as f:
            f.write('\n'.join(accList))
        resdf.to_csv("query/queryResults.csv")

        if datetime.datetime.now().hour >= 16:
            print("Time's up")
            break


if __name__ == "__main__":
    main()