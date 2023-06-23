import numpy as np
import pandas as pd
import statistics as sts

targets = ['16S', '23S', 'atpD', 'groL', 'tuf', 'rpoB', 'concat']
labels = ['16S', '23S', 'atpD', 'hsp65', 'EF-Tu', 'rpoB', 'concatenated sequence']

def snpList(x):
    slist = []
    for ind, row in x.iterrows():
        slist = [*slist, *row.values.tolist()[(2+ind):]]
    return slist


def main():
    data = {target:snpList(pd.read_csv(f"{target}SNP")) for target in targets}
    snpStats = pd.DataFrame(columns=['target', 'avgSNP', 'medianSNP', 'minSNP', 'maxSNP', 'sdSNP'])
    for target in targets:
        temp = data.get(target)
        snpStats.loc[-1] = [target, sts.mean(temp), sts.median(temp), min(temp), max(temp), sts.stdev(temp)]
        snpStats.index = snpStats.index+1
        snpStats = snpStats.sort_index()
    print(snpStats)
    snpStats.to_csv("snpStats.csv")

if __name__ == "__main__":
    main()