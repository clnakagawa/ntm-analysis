import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



def baseMap(c):
    bmap = {'a': 1, 't': 2, 'c': 3, 'g': 4, '-': 0}
    return bmap.get(c)

# take list of nts at pt and count paired diffs
def getDiffs(bases):
    ct = 0
    for i in range(len(bases)):
        for j in range(i,len(bases)):
            if bases[i] != bases[j]:
                ct += 1
    return ct

def seqPlot(seqs, names, target, query=False):
    vals = [[baseMap(seq[i]) for seq in seqs] for i in range(len(seqs[0]))]
    df = pd.DataFrame(vals, columns=names)
    cmap = sns.color_palette("Pastel2", 5)
    ax = sns.heatmap(df, cmap=cmap, vmin=-0.5, vmax=4.5)
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([0,1,2,3,4])
    colorbar.set_ticklabels(['Gap','A','T','C','G'])
    ax.get_figure().savefig(f"plots/{target}basemap")
    plt.clf()

def indelPlot(seqs, target, query=False):
    vals = [[seq[i] for seq in seqs] for i in range(len(seqs[0]))]
    vals = [len(seqs)-pos.count('-') for pos in vals]
    fig = plt.figure()
    plt.plot(vals)
    plt.savefig(f"plots/{target}indelplot")

def diffPlot(seqs, target, query=False):
    diffs = []
    l = min([len(seq) for seq in seqs])
    for i in range(l):
        diffs.append(getDiffs([seq[i] for seq in seqs]))
    fig = plt.figure()
    plt.bar(range(l), diffs)
    plt.savefig(f"plots/{target}_diffplot")
    plt.clf()

def gapPlot(names, target, query=False):
    spNum = {names[i]:i for i in range(len(names))}
    if query:
        df = pd.read_csv(f"queryData/{names[0]}/{target}_gaps.csv")
    else:
        df = pd.read_csv(f"gapData/{target}_gaps.csv")
    fig = plt.figure()
    for index, row in df.iterrows():
        y = spNum.get(row['sp'])
        plt.plot((row['start'], row['end']), (y, y), color='black')
    plt.ylabel('sequences')
    plt.xlabel('sequence position')
    plt.tick_params(left = False, labelleft = False)
    if query:
        plt.savefig(f"queryPlots/{names[0]}/{target}_gapplot")
    else:
        plt.savefig(f"plots/{target}_gapplot")
    plt.clf()

def gapData(seqs, names, target, query=False):
    gapDF = pd.DataFrame(columns=['sp','start','end','len'])
    t = len(seqs)
    for i in range(t):
        sp = names[i]
        seq = seqs[i]
        seql = len(seq)
        start = 0
        end = 0
        for p in range(seql):
            if seq[p] == '-':
                if start:
                    end += 1
                else:
                    start = p+1
                    end = p+2
            else:
                if start:
                    gapDF.loc[len(gapDF)] = [sp, start, end, end-start]
                    start = 0
    if query:
        gapDF.to_csv(f"queryData/{names[0]}/{target}_gaps.csv", index=False)
    else:
        gapDF.to_csv(f"gapData/{target}_gaps.csv", index=False)


def main(targets):
    for target in targets:
        with open(f'targetAligns/{target}_align', 'r') as f:
            seqs = f.read().split("\n>")
        names = [seq.split('\n')[0].replace('>','') for seq in seqs]
        seqs = [''.join(seq.split("\n")[1:]) for seq in seqs]
        diffPlot(seqs, target)
        seqPlot(seqs, names, target)
        indelPlot(seqs, target)
        gapData(seqs, names, target)
        gapPlot(seqs,names,target)



if __name__ == "__main__":
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf', 'ITS', 'concat']
    main(targets)