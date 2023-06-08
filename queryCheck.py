import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import snpDists, seqID, concatenate, alignCheck, os

def main():
    samples = ['test1']
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
    print("gettings SNP and % ID data")
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf', 'concat']
    snpDists.query(targets, samples)
    seqID.query(targets, samples)

    # compile results into one csv
    for sample in samples:
        snps = pd.DataFrame({f"{target}SNP":pd.read_csv(f"queryData/{sample}/{target}SNP.csv").transpose()[0] for target in targets})
        pids = pd.DataFrame({f"{target}PID":pd.read_csv(f"queryData/{sample}/{target}PID.csv").transpose()[0] for target in targets})
        snps.to_csv(f"queryData/{sample}/allSNP.csv")
        pids.to_csv(f"queryData/{sample}/allPID.csv")
        all = snps.join(pids)
        all = all.sort_values(by='concatSNP')
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
        all.to_csv(f"queryData/{sample}/resAll.csv")

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
            plt.clf()



if __name__ == "__main__":
    main()