import numpy as np
import os
import pandas as pd
import re

def main():
    ftpList = pd.read_csv("summaryDataMore.csv")
    subPath = os.path.join(os.getcwd(), "cdsData")
    protList = []
    p = re.compile("\[protein=(.*?)\]")
    for sp in ftpList['sp']:
        fname = f"{sp}_cds.fna"
        with open(os.path.join(subPath, fname), 'r') as f:
            protList += p.findall(f.read())
    protList = list(dict.fromkeys(protList))

    for sp in ftpList['sp']:
        fname = f"{sp}_cds.fna"
        with open(os.path.join(subPath, fname), 'r') as f:
            seqs = f.read()
            for prot in protList:
                if f"[protein={prot}]" not in seqs:
                    protList.remove(prot)

    with open("coreProts.txt", 'w') as f:
        f.write('\n'.join(protList))

if __name__ == "__main__":
    main()