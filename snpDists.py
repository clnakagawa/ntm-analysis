import pandas as pd
import numpy as np

def clean(s):
    cln = ''
    for c in s.upper():
        if c not in ['A', 'T', 'C', 'G']:
            cln += '-'
        else:
            cln += c
    return cln

# calc snp distance btwn 2 plain seqs
def distance(s1, s2):
    cln1 = clean(s1)
    cln2 = clean(s2)
    diff = 0
    for i in range(len(s1)):
        if cln1[i] != cln2[i] and cln2 != '-' and cln1 != '-':
            diff += 1
    return diff

# make pandas df w snp dists labeled by sp
# input is plain fasta
def distMat(fa, target):
    print(target)
    fas = [f.split("\n") for f in fa.split("\n>")]
    names = [f[0].replace(">",'') for f in fas]
    seqs = [''.join(f[1:]) for f in fas]
    mat = [[distance(s1, s2) for s2 in seqs] for s1 in seqs]
    mat = pd.DataFrame(mat, columns=names)
    mat['names'] = names
    mat = mat.set_index('names', drop=True).rename_axis(None)
    print(mat.head())
    mat.to_csv(f"targetSNPs/{target}SNP.csv")

def query(targets, qnames):
    for qname in qnames:
        for target in targets:
            with open(f"queryAligns/{target}_align") as f:
                data = f.read()
            fas = [f.split("\n") for f in data.split("\n>")]
            names = [f[0].replace(">", '') for f in fas]
            qind = names.index(qname)
            seqs = [''.join(f[1:]) for f in fas]
            mat = [[distance(seqs[qind], s) for s in seqs]]
            df = pd.DataFrame(mat, columns=names)
            df = df.drop(qnames, axis=1)
            df.to_csv(f"queryData/{qname}/{target}SNP.csv", index=False)

def main(targets):
    for target in targets:
        with open(f"targetAligns/{target}_align") as f:
            data = f.read()
        distMat(data, target)

if __name__ == "__main__":
    # one target for shorter testing time
    #targets = ['atpD']
    # all targets for updating
    targets = ['atpD', 'groL', 'rpoB', 'tuf', '16S', '23S', 'concat']
    main(targets)