import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



def snpList(x):
    slist = []
    ct = 0
    for ind, row in x.iterrows():
        slist = [*slist, *row.values.tolist()[(2+ct):]]
        ct += 1
    return slist


def main(targets):
    snpdata = {}
    for target in targets:
        print(target)
        fname = f"targetSNPs/{target}SNP.csv"
        snpdata[target] = pd.read_csv(fname, index_col=[0])
        print(snpdata[target].head())
        plt.figure(figsize=(100, 100))
        ax = sns.heatmap(snpdata[target],
                         cmap=sns.color_palette('Reds_r', n_colors=200),
                         square=True,
                         linewidths=0.5,
                         annot=True,
                         fmt='g',
                         annot_kws={"fontsize":8})
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=40)
        plt.savefig(f"plots/{target}SNPmap")
        plt.clf()
    plotData = []
    for target in targets:
        plotData.append(snpList(snpdata.get(target)))
    plt.figure(figsize=(10,7))
    plt.xlabel('Target')
    plt.ylabel('SNP Distance')
    plt.boxplot(plotData, positions=range(1,len(targets)+1), labels=targets)
    plt.savefig("plots/SNPBoxplots")
    plt.show()

if __name__ == "__main__":
    targets = ['16S', '23S', 'rpoB', 'groL', 'atpD', 'tuf', 'concat']
    main(targets)