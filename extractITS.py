import pandas as pd
import numpy as np
import os

def main(makeDB=False):
    pd.set_option('display.max_rows', 200)
    complete = ['Complete Genome', 'Chromosome']
    ftpList = pd.read_csv("summaryData.csv")
    print(ftpList.loc[~ftpList['assembly_level'].isin(complete)]['sp'])
    ftpList = ftpList.loc[ftpList['assembly_level'].isin(complete)]
    print(len(ftpList.index))
    if makeDB:
        for sp in ftpList['sp']:
            fname = f"{sp}_genome.fna"
            print(sp)
            os.system(f"makeblastdb -in genomes/{fname} -dbtype nucl -out blastDBs/{sp}")
            os.system(f"blastn -query refSeqs/ITS.txt -db blastDBs/{sp} -evalue 1e-6 -out itsSeqs/{sp}_its.txt -outfmt \"7 sseq\"")
    # combine all seqs
    allseqs = ''
    for sp in ftpList['sp']:
        with open(f"itsSeqs/{sp}_its.txt", 'r') as f:
            lines = f.read().split('\n')[:-1]
        seqs = [line for line in lines if "#" not in line]
        if seqs:
            seqs.sort(key=len, reverse=True)
            allseqs += f">{sp}\n{seqs[0]}\n"
    with open("targetSeqs/ITS_allseqs.txt", 'w') as f:
        f.write(allseqs)


if __name__ == '__main__':
    main()

