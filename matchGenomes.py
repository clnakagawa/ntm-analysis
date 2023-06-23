import pandas as pd
import numpy as np
import os, re

completes = ['complete', 'Complete']
incompletes = ['contig', 'Contig', 'contig_', 'Contig_', 'Scaffold', 'scaffold', 'Scaffold_', 'scaffold_', 'NODE_']

def getNum(text, pattern):
    if pattern in text:
        part = text.split(pattern)[1]
        match = int(re.findall('\d+', part)[0])
        return match
    return 999

def reorder(x):
    seqs = {seq.split('\n')[0]:''.join(seq.split("\n")[1:]) for seq in x.split("\n>")}
    type = [lab for lab in incompletes if lab in x][0]
    numSeqs = {getNum(key, type):seqs.get(key) for key in seqs.keys()}
    wholeSeq = ''
    for key in sorted(numSeqs.keys()):
        wholeSeq += numSeqs.get(key)
    return wholeSeq

def order(x):
    seqs = [''.join(seq.split("\n")[1:]) for seq in x.split("\n>")]
    return ''.join(seqs)

def main():
    # test with one

    # import all genomes
    ftpList = pd.read_csv("summaryData.csv")
    subPath = os.path.join(os.getcwd(), "genomes")
    toFix = {}
    toWrite = {}
    toOrder = {}
    for sp in ftpList['sp']:
        fname = f"{sp}_genome.fna"
        with open(os.path.join(subPath, fname), 'r') as f:
            content = f.read()
        if any(x in content for x in completes):
            toWrite[sp] = content
        elif any(x in content for x in incompletes):
            toFix[sp] = content
        else:
            toOrder[sp] = content

    print(len(toWrite))
    print(len(toFix))
    print(len(toOrder))

    subPath = os.path.join(os.getcwd(), "genomeSeqs")
    # write reordered seqs
    for key, val in toFix.items():
        fname = f"{key}_genome.txt"
        with open(os.path.join(subPath, fname), 'w') as f:
            f.write(reorder(val))

    # write unordered seqs
    for key, val in toOrder.items():
        fname = f"{key}_genome.txt"
        with open(os.path.join(subPath, fname), 'w') as f:
            f.write(order(val))

    # write complete seqs
    for key, val in toWrite.items():
        fname = f"{key}_genome.txt"
        with open(os.path.join(subPath, fname), 'w') as f:
            f.write(order(val))

    # order contigs correctly

    # rewirte genomes to new folder

    pass

if __name__ == "__main__":
    main()