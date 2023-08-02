import pandas as pd
import numpy as np
import os
import Bio.Align as Align

def removeNum(x):
    return '_'.join(x.split("_")[:-1])

def testClean(targets, name):
    for target in targets:
        if os.path.exists(f"mmRef/{name}/{target}_all.txt"):
            pd.set_option("display.max_rows", 300)
            with open(f"mmRef/{name}/{target}_all.txt", 'r') as f:
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
            with open(f"mmRef/{name}/{target}.txt", 'w') as f:
                f.write(data)

def getSeq(cds, sp, target):
    genes = cds.split("\n>")
    targets = [gene for gene in genes if target in gene]
    targets = [">" + sp + "_" + str(targets.index(gene)+1) + "\n" + "\n".join(gene.split('\n')[1:]) for gene in targets]
    if len(targets) == 0:
        return ""
    seq = "\n".join(targets) + "\n"
    return seq

def test(targets, name, type="cds"):
    for target, id in targets:
        if os.path.exists(f"refFiles/{type}Data/{name}_{type}.fna"):
            subPath = os.path.join(os.getcwd(), f"refFiles/{type}Data")
            fname = f"{name}_{type}.fna"
            with open(os.path.join(subPath, fname), 'r') as f:
                spSeq = getSeq(f.read(), name, target)
            if spSeq == "":
                print(f"Target not found in {name}")
            else:
                print(f"{spSeq.count('>')} targets found in {name}")
            with open(f"mmRef/{name}/{id}_all.txt", 'w') as f:
                f.write(spSeq)

def main():
    # pull target sequences from files
    if not os.path.exists("mmRef"):
        os.mkdir("mmRef")
    tardf = pd.read_csv("targets.csv")
    genetargets = [tuple(r) for r in tardf.loc[tardf['type']=='cds'][['full','short']].to_numpy()]
    print(genetargets)
    rnatargets = [tuple(r) for r in tardf.loc[tardf['type']=='rna'][['full','short']].to_numpy()]
    print(rnatargets)
    targets = tardf['short'].to_list()
    print(targets)
    data = pd.read_csv("summaryData.csv")
    for name in data['sp']:
        if not os.path.exists(f"mmRef/{name}"):
            os.mkdir(f"mmRef/{name}")
        test(genetargets, name)
        test(rnatargets, name, type='rna')
        testClean(targets, name)


if __name__ == "__main__":
    main()
