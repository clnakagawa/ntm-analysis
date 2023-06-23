import requests
import gzip
import os
from pathlib import Path
from requests.exceptions import ConnectTimeout
import numpy as np
import pandas as pd




def getFile(url, spName, out=None, type="genome", use_acc=False):
    spName = spName.replace("[","")
    spName = spName.replace("]","")
    if type != "genome":
        suf = type+"_from_"
    else:
        suf = ""
    ftpPath = url
    acc = ftpPath.split('/')[-1]
    ftpPath = ftpPath + "/" + acc + f"_{suf}genomic.fna.gz"
    if out is None:
        subPath = os.path.join(os.getcwd(), "cdsData")
    else:
        subPath = os.path.join(os.getcwd(), f"{out}/{type}Data")
    filePath = os.path.join(subPath, f"{acc}_{suf}genomic.fna.gz")
    with open(filePath, 'wb') as f:
        for i in range(10):
            try:
                print(ftpPath)
                r = requests.get(ftpPath, timeout=10)
                break
            except:
                print("Timeout, trying again")
        f.write(r.content)
    with gzip.open(filePath, 'rb') as f:
        content = f.read()
    os.remove(filePath)
    if use_acc:
        spName = '_'.join(acc.split('_')[0:2])
    with open(os.path.join(subPath, f"{spName}_{type}.fna"), 'wb') as f:
        f.write(content)
    print(f"File {spName}_{type}.fna written")

def main(file=None, output=None, type=None, acc=False):
    if type is None:
        type = "genome"
    if not os.path.exists(f"{output}/{type}Data"):
        os.makedirs(f"{output}/{type}Data")
    if file is None:
        ftpData = pd.read_csv("summaryData.csv")
    else:
        ftpData = pd.read_csv(file)
    pd.set_option('display.max_columns', 10)
    print(ftpData.head(5))
    ftpData.apply(lambda x: getFile(x['ftp_path'], x['sp'], out=output, type=type, use_acc=acc), axis=1)

def test():
    ftpData = pd.read_csv("testFiles/summaryData.csv")
    pd.set_option('display.max_columns', 10)
    print(ftpData.head(5))
    ftpData.apply(lambda x: getFile(x['ftp_path'], x['# assembly_accession'], out="testFiles"), axis=1)

if __name__ == "__main__":
    main(type="cds")