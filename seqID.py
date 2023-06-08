import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def pid(seq1, seq2):
    ct = 0
    bct = 0
    l = min(len(seq1), len(seq2))
    for i in range(l):
        if seq1[i] == seq2[i]:
            if seq1[i] == '-':
                l = l-1
            else:
                ct += 1
    return ct / l

def query(targets, qnames):
    for qname in qnames:
        for target in targets:
            with open(f"queryAligns/{target}_align") as f:
                data = f.read()
            fas = [f.split("\n") for f in data.split("\n>")]
            names = [f[0].replace(">", '') for f in fas]
            qind = names.index(qname)
            seqs = [''.join(f[1:]) for f in fas]
            mat = [[pid(seqs[qind], s) for s in seqs]]
            df = pd.DataFrame(mat, columns=names)
            df = df.drop(qname, axis=1)
            df.to_csv(f"queryData/{qname}/{target}PID.csv", index=False)

def main(targets):
    for target in targets:
        with open(f"targetAligns/{target}_align") as f:
            seqs = f.read().split("\n>")
        names = [seq.split("\n")[0].replace(">",'') for seq in seqs]
        seqs = [''.join(seq.split("\n")[1:]) for seq in seqs]
        df = pd.DataFrame(columns=['seq1'] + names)

        for i in range(len(seqs)):
            vals = []
            for j in range(len(seqs)):
                vals.append(pid(seqs[i], seqs[j]))
            df.loc[-1] = [names[i]] + vals
            df.index = df.index + 1
            df = df.sort_index()
        print(df.head())
        df = df.set_index('seq1', drop=True).rename_axis(None)
        df.to_csv(f"idTables/{target}_seqid.csv", index=False)
        hmap = sns.heatmap(df)
        fig = hmap.get_figure()
        fig.savefig(f"plots/{target}PIDmap")
        plt.clf()


if __name__ == "__main__":
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf', 'ITS', 'concat']
    main(targets)