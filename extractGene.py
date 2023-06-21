import numpy as np
import pandas as pd
import os
import sys



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


def main(targets, file=None, dir=""):
    if dir != "":
        dir = dir+"/"
    if not os.path.exists(f"{dir}targetSeqs"):
        os.makedirs(f"{dir}targetSeqs")
    for target, id in targets:
        if file == None:
            ftpList = pd.read_csv("summaryData.csv")
        else:
            ftpList = pd.read_csv(file)
        subPath = os.path.join(os.getcwd(), f"{dir}cdsData")
        allSeqs = ""
        for sp in ftpList['sp']:
            fname = f"{sp}_cds.fna"
            with open(os.path.join(subPath, fname), 'r') as f:
                spSeq = getSeq(f.read(), sp, target)
                if spSeq == "":
                    print(f"Target not found in {sp}")
                else:
                    print(f"{spSeq.count('>')} targets found in {sp}")
                    allSeqs += spSeq
        with open(f"{dir}targetSeqs/{id}_allseqs.txt", 'w') as f:
            f.write(allSeqs)

def test(targets, acc):
    for target, id in targets:
        subPath = os.path.join(os.getcwd(), "testFiles/cdsData")

        fname = f"{acc}_cds.fna"
        with open(os.path.join(subPath, fname), 'r') as f:
            spSeq = getSeq(f.read(), acc, target)
        if spSeq == "":
            print(f"Target not found in {acc}")
        else:
            print(f"{spSeq.count('>')} targets found in {acc}")
        with open(f"query/{acc}/{id}_all.txt", 'w') as f:
            f.write(spSeq)

if __name__ == "__main__":
    targets = [('chaperonin GroEL', 'groL'),
               ('ATP synthase subunit beta', 'atpD'),
               ('RNA polymerase subunit beta', 'rpoB'),
               ('elongation factor Tu', 'tuf')]
    main(targets)