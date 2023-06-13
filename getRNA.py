import requests
import gzip
import os
from requests.exceptions import ConnectTimeout
from pathlib import Path
import numpy as np
import pandas as pd


def getFile(url, spName, test=False):
    ftpPath = url
    acc = ftpPath.split('/')[-1]
    ftpPath = ftpPath + "/" + acc + "_rna_from_genomic.fna.gz"
    if test:
        subPath = os.path.join(os.getcwd(), "testFiles/rnaData")
    else:
        subPath = os.path.join(os.getcwd(), "rnaData")
    filePath = os.path.join(subPath, f"{acc}_rna_from_genomic.fna.gz")
    with open(filePath, 'wb') as f:
        for i in range(10):
            try:
                r = requests.get(ftpPath, timeout=10)
                break
            except:
                print("Timeout, trying again")
        f.write(r.content)
    with gzip.open(filePath, 'rb') as f:
        content = f.read()
    os.remove(filePath)
    with open(os.path.join(subPath, f"{spName}_rna.fna"), 'wb') as f:
        f.write(content)
    print(f"File {spName}_rna.fna written")

def main():
    ftpData = pd.read_csv("summaryData.csv")
    pd.set_option('display.max_columns', 10)
    ftpData.apply(lambda x: getFile(x['ftp_path'], x['sp']), axis=1)

def test():
    ftpData = pd.read_csv("testFiles/summaryData.csv")
    pd.set_option('display.max_columns', 10)
    ftpData.apply(lambda x: getFile(x['ftp_path'], x['# assembly_accession'], test=True), axis=1)


if __name__ == "__main__":
    main()