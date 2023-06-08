import numpy as np
import pandas as pd
import os
import sys

targets = [('16S','16S'),('23S','23S')]

# take full fasta input and return seqs with labels
# multiple copies labeled according to input order
def getSeq(cds, sp, target):
    genes = cds.split("\n>")
    targets = [gene for gene in genes if target in gene]
    targets = [">" + sp + "_" + str(targets.index(gene)+1) + "\n" + "\n".join(gene.split('\n')[1:]) for gene in targets]
    if len(targets) == 0:
        return ""
    seq = "\n".join(targets) + "\n"
    return seq


def main(targets):
    for target, id in targets:
        ftpList = pd.read_csv("summaryData.csv")
        subPath = os.path.join(os.getcwd(), "rnaData")
        allSeqs = ""
        for sp in ftpList['sp']:
            fname = f"{sp}_rna.fna"
            with open(os.path.join(subPath, fname), 'r') as f:
                spSeq = getSeq(f.read(), sp, target)
                if spSeq == "":
                    print(f"Target not found in {sp}")
                else:
                    print(f"{spSeq.count('>')} targets found in {sp}")
                    allSeqs += spSeq
        with open(f"targetSeqs/{id}_allseqs.txt", 'w') as f:
            f.write(allSeqs)

if __name__ == "__main__":
    targets = [('16S', '16S'), ('23S', '23S')]
    main(targets)