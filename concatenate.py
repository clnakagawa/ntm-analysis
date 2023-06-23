import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def query(targets, qnames):
    df = pd.read_csv("summaryData.csv")
    names = df['sp']
    data = {qname:'' for qname in qnames}
    for name in names:
        data[name] = ""
    for target in targets:
        fname = f"queryAligns/{target}_align"
        with open(fname, 'r') as f:
            seqs = f.read().split("\n>")
        for seq in seqs:
            name = seq.split("\n")[0].replace(">","")
            seq = ''.join(seq.split("\n")[1:])
            data.update({name:data.get(name)+seq})
    finalSeqs = [f">{key}\n{value}\n" for key, value in data.items()]
    seqlens = [len(''.join(seq.split("\n")[1:])) for seq in finalSeqs]
    finalSeqs = [seq for seq in finalSeqs if len(''.join(seq.split("\n")[1:])) == max(seqlens)]
    content = ''.join(finalSeqs)
    with open("queryAligns/concat_align", 'w') as f:
        f.write(content)


def main(targets):
    df = pd.read_csv("summaryData.csv")
    names = df['sp']
    data = {}
    for name in names:
        data[name] = ""
    for target in targets:
        fname = f"targetAligns/{target}_align"
        with open(fname, 'r') as f:
            seqs = f.read().split("\n>")
        for seq in seqs:
            name = seq.split("\n")[0].replace(">","")
            seq = ''.join(seq.split("\n")[1:])
            data.update({name:data.get(name)+seq})
    finalSeqs = [f">{key}\n{value}\n" for key, value in data.items()]
    seqlens = [len(''.join(seq.split("\n")[1:])) for seq in finalSeqs]
    finalSeqs = [seq for seq in finalSeqs if len(''.join(seq.split("\n")[1:])) == max(seqlens)]
    seqlens = [l for l in seqlens if l == max(seqlens)]
    content = ''.join(finalSeqs)
    with open("targetAligns/concat_align", 'w') as f:
        f.write(content)


if __name__ == "__main__":
    targets = ['atpD', 'groL', 'rpoB', 'tuf', '16S', '23S']
    main(targets)