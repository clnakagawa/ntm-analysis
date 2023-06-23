import pandas as pd
from Bio import Align
import cleanFTP, getGenes, random, seqClean, os, snpDists, datetime, re
import extractGene
from functools import reduce

def test(sample, source="refFiles"):
    start = datetime.datetime.now()
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf']
    if not os.path.exists(f"query/{sample}"):
        os.mkdir(f"query/{sample}")
        # pull sequences and send to query
        genetargets = [('chaperonin GroEL', 'groL'),
                   ('ATP synthase subunit beta', 'atpD'),
                   ('RNA polymerase subunit beta', 'rpoB'),
                   ('elongation factor Tu', 'tuf')]


        extractGene.test(genetargets, sample)
        rnatargets = [('16S', '16S'), ('23S', '23S')]
        extractGene.test(rnatargets, sample, type='rna')
        seqClean.test(targets, sample)
    dfs = []
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2

    for target in targets:
        with open(f"query/{sample}/{target}.txt") as f:
            qseq = ''.join(f.read().split("\n")[1:])
        with open(f"{source}/cleanedSeqs/{target}.txt", 'r') as f:
            refs = f.read().split("\n>")
        seqs = [''.join(seq.split("\n")[1:]) for seq in refs]
        sps = [seq.split("\n")[0].replace(">", "") for seq in refs]
        alns = [aligner.align(qseq, seq)[0] for seq in seqs]

        snps = snpDists.distanceList([aln[0] for aln in alns],
                                     [aln[1] for aln in alns])
        gaps = [ctgaps(aln[0]) + ctgaps(aln[1]) for aln in alns]
        df = pd.DataFrame({'sp': sps, 'snp': snps, 'gap': gaps}).sort_values('snp')
        df = df.set_index(['sp'])
        df.to_csv(f"query/{sample}/{target}SNP{source}.csv")
        dfs.append(df)
    total = reduce(lambda x, y: x.add(y), dfs)
    total = total.sort_values('gap')
    total = total.sort_values('snp')
    total.to_csv(f"query/{sample}/totalSNP{source}.csv")
    print(total.head())
    print(sample,
          total.index.values[0],
          total['snp'].values[0],
          total['gap'].values[0]
          )
    print(f"Time is {datetime.datetime.now() - start}")



def ctgaps(seq):
    ct = 0
    current = False
    for c in seq:
        if c == '-':
            if not current:
                ct += 1
                current = True
        else:
            if current:
                current = False
    return ct

def main(source="refFiles"):
    pd.set_option('display.max_columns', 15)
    pd.set_option('display.max_rows', 400)
    # cleanFTP.test()
    # getGenes.main(file="testFiles/summaryData.csv", output="testFiles", type="cds", acc=True)
    # getGenes.main(file="testFiles/summaryData.csv", output="testFiles", type="rna", acc=True)
    # getGenes.main(file="testFiles/summaryData.csv", output="testFiles", acc=True)

    # check if all targets present
    # get accession numbers
    # OLD CODE from getting the list of ones w targets
    sumData = pd.read_csv("testFiles/summaryData.csv")
    # accnums = sumData['# assembly_accession']
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
        with open(f"query/toTest{source}.txt", 'r') as f:
            accList = f.read().split("\n")
    except:
        accList = list(df['accs'])

    try:
        resdf = pd.read_csv(f'query/queryResults{source}.csv', index_col=0)
    except:
        resdf = pd.DataFrame(columns=['acc',
                                      'name',
                                      'id',
                                      'realDist',
                                      'idDist',
                                      'realGaps',
                                      'idGaps'])
    try:
        sumdf = pd.read_csv(f'query/querySummary{source}.csv', index_col=0)
    except:
        sumdf = pd.DataFrame(columns=['sp','accuracy','avgDist','frequency'])
    try:
        errordf = pd.read_csv(f'query/queryErrors{source}.csv', index_col=0)
    except:
        errordf = pd.DataFrame(columns=['actual', 'result', 'frequency'])



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
        extractGene.test(rnatargets, acc, type='rna')
        seqClean.test(targets, acc)
        dfs = []
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2
        print(aligner.mismatch_score)

        for target in targets:
            with open(f"query/{acc}/{target}.txt") as f:
                qseq = ''.join(f.read().split("\n")[1:])
            with open(f"{source}/cleanedSeqs/{target}.txt", 'r') as f:
                refs = f.read().split("\n>")
            seqs = [''.join(seq.split("\n")[1:]) for seq in refs]
            sps = [seq.split("\n")[0].replace(">","") for seq in refs]
            alns = [aligner.align(qseq, seq)[0] for seq in seqs]
            snps = snpDists.distanceList([aln[0] for aln in alns],
                                         [aln[1] for aln in alns])
            gaps = [ctgaps(aln[0]) + ctgaps(aln[1]) for aln in alns]
            df = pd.DataFrame({'sp':sps,'snp':snps,'gap':gaps}).sort_values('snp')
            df = df.set_index(['sp'])
            df.to_csv(f"query/{acc}/{target}SNP{source}.csv")
            dfs.append(df)
        total = reduce(lambda x, y: x.add(y), dfs)
        total = total.sort_values('gap')
        total = total.sort_values('snp')
        total.to_csv(f"query/{acc}/totalSNP{source}.csv")
        print(total.head())
        realEntry = sumData.loc[sumData['# assembly_accession'] == acc]
        realSp = realEntry['sp'].values[0]

        # handle case where real species isn't included
        if realSp not in total.index:
            print("actual species not on file")
            realDist = 0
            realGap = 0
        else:
            realDist = total.loc[realEntry['sp'].values[0]]['snp']
            realGap = total.loc[realEntry['sp'].values[0]]['gap']

        print(acc,
              realSp,
              total.index.values[0],
              realDist,
              total['snp'].values[0],
              realGap,
              total['gap'].values[0]
              )
        resdf.loc[len(resdf.index)] = [acc,
                                       realSp,
                                       total.index.values[0],
                                       realDist,
                                       total['snp'].values[0],
                                       realGap,
                                       total['gap'].values[0]
                                       ]
        sp = cleanFTP.getSp(realEntry['organism_name'].values[0])
        # update summary df
        current = sumdf.loc[sumdf['sp'] == sp]
        if not current.empty:
            avg = current['accuracy'].values[0]
            ct = current['frequency'].values[0]
            dist = current['avgDist'].values[0]
            sumdf.loc[sumdf['sp'] == sp] = [sp,
                                             (avg*ct+int(sp==total.index.values[0]))/(ct+1),
                                             (dist*ct+realDist)/(ct+1),
                                             ct+1]
        else:
            sumdf.loc[len(sumdf.index)] = [sp,
                                           int(sp == total.index.values[0]),
                                           realDist,
                                           1]
        print(sumdf.head())

        # update error df
        errordf.loc[len(errordf.index)] = [sp, total.index.values[0], 1]
        errordf = errordf.groupby(['actual', 'result'])['frequency'].sum().reset_index()
        print(errordf.head())
        print(f"Time is {datetime.datetime.now()-start}")
        with open(f"query/toTest{source}.txt", 'w') as f:
            f.write('\n'.join(accList))
        # update all results csv
        resdf.to_csv(f"query/queryResults{source}.csv")
        # update summary csv
        sumdf.to_csv(f'query/querySummary{source}.csv')
        # update error csv
        errordf.to_csv(f'query/queryErrors{source}.csv')
        if datetime.datetime.now().hour >= 16:
            print("Time's up")
            break


if __name__ == "__main__":
    main()
    main(source="lpsnFiles")