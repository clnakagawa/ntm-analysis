import pandas as pd
import numpy as np
from itertools import combinations

def getT(seq):
    t = 0
    for b in seq:
        if b in ['a', 't']:
            t += 2
        else:
            t += 4
    return t

def getSubseq(reg):
    length = len(reg) + 1
    return [reg[x:y] for x, y, in combinations(range(length), r=2) if y > x + 19]

# search conserved region for primers
def searchRegion(regs):
    all = []
    for reg in regs:
        new = getSubseq(reg[0])
        all += [(seq, reg[1]) for seq in new]
    return all

# search aligned sequences to find conserved region
def findRegion(seqs):
    print(len(seqs[0]))
    seqs = [[seq[i] for seq in seqs] for i in range(len(seqs[0]))]
    t = ''
    regs = []
    for b in seqs:
        sb = set(b)
        if len(sb) < 2:
            t += list(sb)[0]
        else:
            if len(t)>19:
                regs.append((t, seqs.index(b)))
            t = ''
    print(regs)
    return regs

def primers(seqs):
    #data = pd.DataFrame(columns=['seq', 'len', 'start', 'end','Tm'])
    regions = findRegion(seqs)
    data = [[seq[0], len(seq[0]), seq[1], seq[1]+len(seq[0]), getT(seq[0])]for seq in searchRegion(regions)]
    data = pd.DataFrame(data, columns=['seq', 'len', 'start', 'end', 'Tm'])
    print(data.head())
    return data


def main(targets):
    for target in targets:
        print(target)
        with open(f"targetAligns/{target}_align", 'r') as f:
            seqs = f.read().split("\n>")
            seqs = [''.join(seq.split("\n")[1:]) for seq in seqs]
            res = primers(seqs)
            res.to_csv(f"primerData/{target}Primers.csv")

if __name__ == "__main__":
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf']
    main(targets)