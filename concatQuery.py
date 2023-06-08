import pandas as pd
import numpy as np

targets = ['atpD', 'groL', 'rpoB', 'tuf', '16S', '23S']


def main():
    df = pd.read_csv("summaryData.csv")
    names = list(df['sp'])
    names.append('test1')
    data = {}
    for name in names:
        data[name] = ""
    for target in targets:
        fname = f"{target}test_align"
        with open(fname, 'r') as f:
            seqs = f.read().split("\n>")
        for seq in seqs:
            name = seq.split("\n")[0].replace(">","")
            seq = ''.join(seq.split("\n")[1:])
            data.update({name:data.get(name)+seq})
    finalSeqs = [f">{key}\n{value}\n" for key, value in data.items()]
    finalSeqs = [seq for seq in finalSeqs if len(seq) > 17500]
    content = ''.join(finalSeqs)
    with open("concatTestAlign.txt", 'w') as f:
        f.write(content)


if __name__ == "__main__":
    main()