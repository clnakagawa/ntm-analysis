from Bio import Align
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

lowLim = 0.95
upLim = 1.25

def test(targets, acc):
    for target in targets:
        pd.set_option("display.max_rows", 300)
        with open(f"query/{acc}/{target}_all.txt", 'r') as f:
            testSeqs = f.read().split(">")[1:]
        seqData = {'name':[x.split("\n")[0] for x in testSeqs],
                   'seq':[''.join(x.split('\n')[1:]) for x in testSeqs]}
        seqData = pd.DataFrame(seqData)

        with open(f"refSeqs/{target}.txt") as f:
            seq1 = ''.join(f.read().split("\n")[1:])

        aligner = Align.PairwiseAligner()
        seqData['score'] = seqData['seq'].apply(lambda x: aligner.score(seq1, x))
        seqData['score'] = pd.to_numeric(seqData['score'])


        seqData = seqData.sort_values("score", ascending=False)
        seqData['sp'] = seqData['name'].apply(lambda x: removeNum(x))
        seqData = seqData.drop_duplicates(subset=['sp'])
        seqData['fa'] = seqData.apply(lambda x: ">" + x['sp'] + "\n" + x['seq'], axis=1)
        data = '\n'.join(seqData['fa'])
        with open(f"query/{acc}/{target}.txt", 'w') as f:
            f.write(data)


def removeNum(x):
    return '_'.join(x.split("_")[:-1])

def main(targets, test=False):
    for target in targets:
        pd.set_option("display.max_rows", 300)
        with open(f"refSeqs/{target}.txt", 'r') as f:
            seq1 = "".join(f.read().split("\n")[1:])
        if test:
            with open(f"query/{target}_all.txt", 'r') as f:
                testSeqs = f.read().split(">")[1:]
        else:
            with open(f"targetSeqs/{target}_allseqs.txt", 'r') as f:
                testSeqs = f.read().split(">")[1:]
        seqData = {'name':[x.split("\n")[0] for x in testSeqs],
                   'seq':[''.join(x.split('\n')[1:]) for x in testSeqs]}
        seqData = pd.DataFrame(seqData)

        if not test:
        # filter short seqs here
            seqData['len'] = seqData['seq'].str.len()
            seqData = seqData.loc[seqData['len'] > len(seq1) * lowLim]
            seqData = seqData.loc[seqData['len'] < len(seq1) * upLim]

        aligner = Align.PairwiseAligner()
        seqData['score'] = seqData['seq'].apply(lambda x: aligner.score(seq1, x))
        seqData['score'] = pd.to_numeric(seqData['score'])


        seqData = seqData.sort_values("score", ascending=False)
        seqData['sp'] = seqData['name'].apply(lambda x: removeNum(x))
        seqData = seqData.drop_duplicates(subset=['sp'])

        # check scores again
        avgScore = np.mean(seqData['score'])
        sdScore = np.std(seqData['score'])
        if not test:
            seqData = seqData.loc[seqData['score'] > avgScore - 3 * sdScore]

        seqData.hist(column='score')
        pltname = target + "_scores.png"
        plt.savefig(f"plots/{pltname}")

        seqData['fa'] = seqData.apply(lambda x: ">" + x['sp'] + "\n" + x['seq'], axis=1)
        with open(f"query/{target}.txt", 'w') as f:
            f.write('\n'.join(seqData['fa']))

if __name__ == "__main__":
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf']
    main(targets)