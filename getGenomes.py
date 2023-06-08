import requests
import gzip
import os
from pathlib import Path
import numpy as np
import pandas as pd


def getFile(url, spName):
    ftpPath = url
    acc = ftpPath.split('/')[-1]
    ftpPath = ftpPath + "/" + acc + "_genomic.fna.gz"
    subPath = os.path.join(os.getcwd(), "genomes")
    filePath = os.path.join(subPath, f"{acc}_genomic.fna.gz")
    with open(filePath, 'wb') as f:
        r = requests.get(ftpPath)
        f.write(r.content)
    with gzip.open(filePath, 'rb') as f:
        content = f.read()
    os.remove(filePath)
    with open(os.path.join(subPath, f"{spName}_genome.fna"), 'wb') as f:
        f.write(content)
    print(f"File {spName}_genome.fna written")

def main():
    ftpData = pd.read_csv("summaryData.csv")
    pd.set_option('display.max_columns', 10)
    print(ftpData.head(5))
    ftpData.apply(lambda x: getFile(x['ftp_path'], x['sp']), axis=1)


if __name__ == "__main__":
    main()