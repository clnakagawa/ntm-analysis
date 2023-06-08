import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



def snpList(x):
    slist = []
    for ind, row in x.iterrows():
        slist = [*slist, *row.values.tolist()[(2+ind):]]
    return slist


def main(targets):
    snpdata = {}
    for target in targets:
        fname = f"targetSNPs/{target}SNP"
        snpdata[target] = pd.read_csv(fname)
        hmap = sns.heatmap(snpdata[target])
        fig = hmap.get_figure()
        fig.savefig(f"plots/{target}SNPmap")
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
    targets = ['16S', '23S', 'rpoB', 'groL', 'atpD', 'tuf', 'ITS', 'concat']
    main(targets)