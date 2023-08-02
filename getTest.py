import pandas as pd
import numpy as np
from Bio import Align
import matplotlib.pyplot as plt
import cleanFTP, getGenes, random, seqClean, os, snpDists, datetime, re
import extractGene
from functools import reduce

def gapDiff(seq1, seq2):
    ct = 0
    for i in range(len(seq1)):
        if (seq1[i] != '-') and (seq2[i] == '-'):
            ct += 1
    return ct 

def getCov(seq1, seq2):
    ct = 0
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            ct += 1
    return ct/(len(seq1)-seq1.count("-"))


def test(sample, source="refFiles"):
    # must change if targets are changed
    targets = pd.read_csv('targets.csv')['short'].to_list()
    print(targets)
    start = datetime.datetime.now()
    if os.path.exists(f"query/{sample}/totalSNP{source}.csv"):
        total = pd.read_csv(f"query/{sample}/totalSNP{source}.csv", index_col=0)
    else:
        if not os.path.exists(f"query/{sample}"):
            os.mkdir(f"query/{sample}")
        dfs = []
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2
        with open(f"query/{sample}/allSeqs.txt", 'r') as f:
            fseqs = f.read().split("\n>")
            qseqs = [''.join(seq.split("\n")[1:]) for seq in fseqs]
        with open(f"query/{sample}/log.txt", 'w') as f:
            f.write("")
        for target in targets:
            if not os.path.exists(f"query/{sample}/{target}"):
                os.mkdir(f"query/{sample}/{target}")
            qseq = qseqs[targets.index(target)]
            with open(f"query/{sample}/log.txt", 'a') as f:
                f.write("\n" + fseqs[targets.index(target)] + "\nFullseqs\n" + qseqs[targets.index(target)])
            with open(f"{source}/cleanedSeqs/{target}.txt", 'r') as f:
                refs = f.read().split("\n>")
            seqs = [''.join(seq.split("\n")[1:]) for seq in refs]
            sps = [seq.split("\n")[0].replace(">", "") for seq in refs]
            alns = [aligner.align(qseq, seq)[0] for seq in seqs]
            for i in range(len(alns)):
                with open(f"query/{sample}/{target}/align{target}_{sps[i]}", 'w') as f:
                    f.write(alns[i].format('fasta'))
            snps = snpDists.distanceList([aln[0] for aln in alns],
                                         [aln[1] for aln in alns])
            gaps = [ctgaps(aln[0]) + ctgaps(aln[1]) for aln in alns]
            cov = [getCov(aln[0],aln[1]) for aln in alns]
            gd = [gapDiff(aln[0],aln[1]) for aln in alns]
            df = pd.DataFrame({'sp': sps, 'snp': snps, 'gap': gaps, 'cov': cov, 'gdiff': gd})
            df['diff'] = df['snp'] + df['gdiff']
            df = df.sort_values("diff")
            df = df.set_index(['sp'])
            df.to_csv(f"query/{sample}/{target}SNP{source}.csv")
            dfs.append(df)
        total = reduce(lambda x, y: x.add(y), dfs)
        total['cov'] = total['cov'] / 6
        total = total.sort_values('gap')
        total = total.sort_values('diff')
        total.to_csv(f"query/{sample}/totalSNP{source}.csv")
    print(total.head())
    print(sample,
          total.index.values[0],
          total['snp'].values[0],
          total['gap'].values[0],
          total['cov'].values[0],
          total['gdiff'].values[0]
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

def ctmatches(seq1, seq2):
    ct = 0
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            ct += 1
    return ct

def main(source="refFiles", redo=False):
    pd.set_option('display.max_columns', 15)
    pd.set_option('display.max_rows', 400)
    
    if redo:
        cleanFTP.test()
        getGenes.main(file="testFiles/summaryData.csv", output="testFiles", type="cds", acc=True)
        getGenes.main(file="testFiles/summaryData.csv", output="testFiles", type="rna", acc=True)
        getGenes.main(file="testFiles/summaryData.csv", output="testFiles", acc=True)

    # check if all targets present
    # get accession numbers
    # OLD CODE from getting the list of ones w targets
    sumData = pd.read_csv("testFiles/summaryData.csv")
    if redo:
        accnums = sumData['# assembly_accession']
        # must change if targets change
        geneNames = ['ATP synthase subunit beta', 'chaperonin GroEL', 'RNA polymerase subunit beta', 'elongation factor Tu']
        rnaNames = ['16S', '23S']
    
        for target in geneNames:
            print(target)
            acctemp = []
            for a in accnums:
                with open(f"testFiles/cdsData/{a}_cds.fna", 'r') as f:
                    if target in f.read():
                        acctemp.append(a)
            accnums = acctemp
        for target in rnaNames:
            print(target)
            acctemp = []
            for a in accnums:
                with open(f"testFiles/rnaData/{a}_rna.fna", 'r') as f:
                    if target in f.read():
                        acctemp.append(a)
            accnums = acctemp
    
        print(len(accnums))
        df = pd.DataFrame({'accs':accnums})
        df.to_csv("testFiles/hasTargetList.csv")

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
        sumgendf = pd.read_csv(f'query/querySummaryGen{source}.csv', index_col=0)
    except:
        sumgendf = pd.DataFrame(columns=['sp','accuracy','avgDist','frequency','top_hits'])

    try:
        errordf = pd.read_csv(f'query/queryErrors{source}.csv', index_col=0)
    except:
        errordf = pd.DataFrame(columns=['actual', 'result', 'frequency'])


    while accList:
        acc = accList.pop()
        start = datetime.datetime.now()
        if os.path.exists(f"query/{acc}/totalSNP{source}.csv"):
            total = pd.read_csv(f"query/{acc}/totalSNP{source}.csv", index_col=0)
        else:
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
            aligner.mismatch_score = -2

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
                matches = [ctmatches(aln[0], aln[1]) for aln in alns]
                for i in range(len(alns)):
                    with open(f"query/{acc}/align{target}_{sps[i]}",'w') as f:
                        f.write(alns[i].format('fasta'))
                df = pd.DataFrame({'sp':sps,'snp':snps,'gap':gaps,'match':matches}).sort_values('snp')
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
        print(realSp)

        # handle case where real species isn't included
        if realSp not in total.index.values:
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

        # update general summary df
        current = sumgendf.loc[sumgendf['sp'] == sp]
        if not current.empty:
            avg = current['accuracy'].values[0]
            ct = current['frequency'].values[0]
            dist = current['avgDist'].values[0]
            hits = current['top_hits'].values[0]
            print(hits)
            if isinstance(hits, float):
                hits = ""
            sumgendf.loc[sumgendf['sp'] == sp] = [sp,
                                               (avg*ct+int(sp in total[total['snp'] < 10].index.values))/(ct+1),
                                               (dist*ct+realDist)/(ct+1),
                                               ct+1,
                                               ';'.join(list(set(list(total[total['snp'] < 10].index.values) + hits.split(";"))))]
        else:
            sumgendf.loc[len(sumgendf.index)] = [sp,
                                                 int(sp in total[total['snp'] < 10].index.values),
                                                 realDist,
                                                 1,
                                                 ';'.join(total[total['snp'] < 10].index.values)]

        # update error df
        errordf.loc[len(errordf.index)] = [sp, total.index.values[0], 1]
        errordf = errordf.groupby(['actual', 'result'])['frequency'].sum().reset_index()
        print(f"Time is {datetime.datetime.now()-start}")
        with open(f"query/toTest{source}.txt", 'w') as f:
            f.write('\n'.join(accList))
        # update all results csv
        resdf.to_csv(f"query/queryResults{source}.csv")
        # update summary csv
        sumdf.to_csv(f'query/querySummary{source}.csv')
        # update summary csv
        sumgendf.to_csv(f'query/querySummaryGen{source}.csv')
        # update error csv
        errordf.to_csv(f'query/queryErrors{source}.csv')

def runTest():
    main()
    main(source="lpsnFiles")

if __name__ == "__main__":
    pass
