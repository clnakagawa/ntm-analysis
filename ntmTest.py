import pandas as pd
import numpy as np
from functools import reduce

targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf', 'concat']

def getHits(target):
    df = pd.read_csv(f"querySNP/{target}TestSNP")
    df = df.sort_values('test1')
    df[target] = df['test1']
    df = df.loc[df['snp-dists 0.8.2'] != 'test1']
    return df[['snp-dists 0.8.2', target]]

def main():
    pd.set_option('display.max_columns', 8)
    dfs = []
    for target in targets:
        dfs.append(getHits(target))
    df = reduce(lambda left,right: pd.merge(left, right, on=['snp-dists 0.8.2'], how='outer'), dfs)
    df = df.sort_values('concat')
    print(df.head(10))

if __name__ == "__main__":
    main()