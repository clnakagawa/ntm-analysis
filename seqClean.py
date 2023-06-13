from Bio import AlignIO, pairwise2, Seq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

lowLim = 0.95
upLim = 1.25



def removeNum(x):
    return '_'.join(x.split("_")[:-1])

def alnScore(s1, s2):
    return pairwise2.align.globalxx(s1, s2, score_only=True)


def main(targets, test=False):
    for target in targets:
        pd.set_option("display.max_rows", 300)
        with open(f"refSeqs/{target}.txt", 'r') as f:
            seq1 = "".join(f.read().split("\n")[1:])
        print(seq1)
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
            print(seqData.loc[seqData['seq'].str.len() < len(seq1) * lowLim])
            seqData = seqData.loc[seqData['len'] > len(seq1) * lowLim]
            print(seqData.loc[seqData['seq'].str.len() > len(seq1) * upLim])
            seqData = seqData.loc[seqData['len'] < len(seq1) * upLim]


        seqData['score'] = seqData['seq'].apply(lambda x: alnScore(seq1, x))
        seqData['score'] = pd.to_numeric(seqData['score'])


        seqData = seqData.sort_values("score", ascending=False)
        seqData['sp'] = seqData['name'].apply(lambda x: removeNum(x))
        seqData = seqData.drop_duplicates(subset=['sp'])

        # check scores again
        avgScore = np.mean(seqData['score'])
        sdScore = np.std(seqData['score'])
        print(avgScore)
        print(sdScore)
        if not test:
            seqData = seqData.loc[seqData['score'] > avgScore - 3 * sdScore]

        seqData.hist(column='score')
        pltname = target + "_scores.png"
        plt.savefig(f"plots/{pltname}")

        print(seqData['name'])
        seqData['fa'] = seqData.apply(lambda x: ">" + x['sp'] + "\n" + x['seq'], axis=1)
        print(len(seqData.index))
        with open(f"query/{target}.txt", 'w') as f:
            f.write('\n'.join(seqData['fa']))

if __name__ == "__main__":
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf']
    main(targets)