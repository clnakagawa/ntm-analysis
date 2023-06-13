import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import snpDists, seqID, concatenate, alignCheck, os

def main():
    with open("testFiles/testList.txt", 'r') as f:
        samples = f.read().split("\n")
    for sample in samples:
        try:
            os.mkdir(f"queryData/{sample}")
            os.mkdir(f"queryPlots/{sample}")
        except OSError as error:
            print(error)
    pd.set_option('display.max_columns', 14)
    # concatenate
    print("concatenating sequences")
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf']
    concatenate.query(targets, samples)

    # run snp-dists and seqId
    # print("gettings SNP and % ID data")
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf', 'concat']
    snpDists.query(targets, samples)
    seqID.query(targets, samples)

    ftpData = pd.read_csv("testFiles/summaryData.csv")

    # open csv for data recording
    try:
        stats = pd.read_csv('testFiles/queryStats.csv', index_col=0)
    except:
        stats = pd.DataFrame(columns=['accuracy','avgSNP','n'])
    print(stats.head())
    # compile results into one csv
    for sample in samples:
        snps = pd.DataFrame({f"{target}SNP":pd.read_csv(f"queryData/{sample}/{target}SNP.csv").transpose()[0] for target in targets})
        pids = pd.DataFrame({f"{target}PID":pd.read_csv(f"queryData/{sample}/{target}PID.csv").transpose()[0] for target in targets})
        snps.to_csv(f"queryData/{sample}/allSNP.csv")
        pids.to_csv(f"queryData/{sample}/allPID.csv")
        all = snps.join(pids)
        all = all.sort_values(by='concatSNP')

        # find top results and correct result
        answer = ftpData.loc[ftpData['# assembly_accession'] == sample]['sp'].values[0]
        print(answer)
        answerDist = all.loc[answer]['concatSNP']
        print(answerDist)
        res = all.index[0]
        print(res)
        resDist = all['concatSNP'].values[0]
        print(resDist)

        if answer in stats.index.values:
            # update
            n = stats.loc[answer]['n']
            accuracy = stats.loc[answer]['accuracy']
            avgSNP = stats.loc[answer]['avgSNP']
            stats.loc[answer] = [(n * accuracy + (answer == res)) / (n + 1),
                                 (n * avgSNP + (answerDist)) / (n + 1),
                                 n+1]
        else:
            stats.loc[answer] = [int(answer == res),
                                 answerDist,
                                 1]


        print(all.head(10))
        for target in targets:
            fig = plt.figure()
            plt.hist(all[f"{target}SNP"])
            plt.xlabel("SNP distance")
            plt.savefig(f"queryPlots/{sample}/{target}SNPDist")
            plt.clf()
            fig = plt.figure()
            plt.hist(all[f"{target}PID"])
            plt.xlabel("% identity")
            plt.savefig(f"queryPlots/{sample}/{target}PIDDist")
            plt.clf()
            plt.close('all')
        all.to_csv(f"queryData/{sample}/resAll.csv")
        with open(f"queryData/{sample}/answer.txt", 'w') as f:
            print(ftpData.loc[ftpData['# assembly_accession'] == sample])
            f.write(answer)
    stats.to_csv("testFiles/queryStats.csv")
    # make plots
    for target in targets:
        with open(f'queryAligns/{target}_align', 'r') as f:
            seqs = [seq.split('\n') for seq in f.read().split("\n>")]
        names = [seq[0].replace('>','') for seq in seqs]
        seqs = [''.join(seq[1:]) for seq in seqs]
        seqs = [seqs[names.index(sample)] for sample in samples]
        names = samples
        for i in range(len(seqs)):
            alignCheck.gapData(seqs[i:i+1], names[i:i+1], target, query=True)
            alignCheck.gapPlot(names[i:i+1], target, query=True)
            df = pd.read_csv(f'queryData/{names[i]}/{target}_gaps.csv')
            fig = plt.figure()
            plt.hist(df['len'])
            plt.xlabel('Gap Length')
            plt.savefig(f"queryPlots/{names[i]}/{target}GapHist")
            plt.close('all')



if __name__ == "__main__":
    main()