import pandas as pd
import numpy as np
from Bio.Align.Applications import MafftCommandline
import os, subprocess
def run():
    # os.chdir("C:\\Users\cyn06\OneDrive - New York State Office of Information Technology Services\Downloads\mafft-7.520-win64-signed\mafft-win")
    # os.startfile('qAlign.bat')
    profile1 = "targetAligns/16S_align"
    newSeqs = "query/16S.txt"
    outfile = "muscleTest.txt"
    os.system(f"muscle --profile -in1 {profile1} -in2 {newSeqs} -out {outfile}")
    musclePath = "C:\\Users\cyn06\OneDrive - New York State Office of Information Technology Services\Downloads\myuscle5.1.win64.exe"

def main():
    run()

if __name__ == "__main__":
    main()