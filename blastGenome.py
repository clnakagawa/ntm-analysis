import numpy as np
import pandas as pd
import os, subprocess, getTest

def main():
    pd.set_option("display.max_columns", 10)
    # find files in genomesTest
    # if matching files in query don't exist, run and create files
    files = [f for f in os.listdir("genomesTest") if os.path.isfile(os.path.join("genomesTest", f))]
    fnames = ['.'.join(f.split(".")[:-1]) for f in files]
    print(fnames)

    for f in fnames:
        if not os.path.exists(f"query/{f}"):
            os.makedirs(f"query/{f}")


    # concatenate seqs
    for i in range(len(files)):
        file = files[i]
        with open(f"genomesTest/{file}", 'r') as f:
            parts = f.read().split("\n>")
        parts = [''.join(part.split("\n")[1:]) for part in parts]
        with open(f"genomesTest/{file}", 'w') as f:
            f.write(">"+fnames[i]+"\n"+''.join(parts))

    # make blastdbs
    for i in range(len(files)):
        if not os.path.exists(f"query/{fnames[i]}/{fnames[i]}.ndb"):
            cmd = [os.environ['BLAST']+"\makeblastdb",
                   '-in',
                   f"genomesTest/{files[i]}",
                   '-dbtype',
                   'nucl',
                   '-out',
                   f'query/{fnames[i]}/{fnames[i]}']
            subprocess.call(cmd, shell=True)

    # query blastdbs
    targets = ['16S', '23S', 'groL', 'atpD', 'rpoB', 'tuf', 'ITS']
    for target in targets:
        for f in fnames:
            cmd = [os.environ['BLAST']+"\\blastn",
                   "-query",
                   f"refSeqs/{target}.txt",
                   "-db",
                   f"query/{f}/{f}",
                   "-evalue",
                   "1e-6",
                   "-out",
                   f"query/{f}/{target}_blast.txt",
                   "-outfmt",
                   "7 sseq length qstart qend sstart send bitscore"]
            subprocess.call(cmd)
            df = pd.read_table(f"query/{f}/{target}_blast.txt", skiprows=5, skipfooter=1, header=None, engine='python')
            df.columns = ['seq', 'len', 'qstart', 'qend', 'sstart', 'send', 'bscore']
            df = df.loc[df['bscore'] > df['len'] / 2]
            df = df.sort_values('qstart')

            # concatenate parts and write to sequence
            # add gaps in where no coverage
            pos = 0
            seq = ""
            for ind, row in df.iterrows():
                seq += '-'*(row['qstart'] - pos)
                seq += row['seq']
                pos = row['qend']

            with open(f"query/{f}/{target}.txt", 'w') as file:
                file.write(f">{f}\n{seq}")

    for f in fnames:
        print(f"Test using refFiles for {f}")
        getTest.test(f, source='refFiles')
        print(f"Test using lpsnFiles for {f}")
        getTest.test(f, source='lpsnFiles')



if __name__ == "__main__":
    main()