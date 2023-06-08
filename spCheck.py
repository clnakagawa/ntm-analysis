import pandas as pd
import numpy as np

# script to look at chimaera vs intracellulare
targets = ['16S', '23S', 'atpD', 'groL', 'atpD', 'tuf', 'concat']

def main():
    data = []
    resDF = pd.DataFrame(columns=['16S_SNP', '16S_%_identity',
                                  '23S_SNP', '23S_%_identity',
                                  'atpD_SNP', 'atpD_%_identity',
                                  'groL_SNP', 'groL_%_identity',
                                  'rpoB_SNP', 'rpoB_%_identity',
                                  'tuf_SNP', 'tuf_%_identity',
                                  'concat_SNP', 'concat_%_identity'])
    for target in targets:
        df = pd.read_csv(f"targetSNPs/{target}SNP")
        print(f"\n{target} SNP distance")
        print(df.loc[df['snp-dists 0.8.2'] == 'Mycobacterium_intracellulare', 'Mycobacterium_chimaera'].values[0])
        data.append(df.loc[df['snp-dists 0.8.2'] == 'Mycobacterium_intracellulare', 'Mycobacterium_chimaera'].values[0])
        df = pd.read_csv(f"idTables/{target}_seqid.csv")
        print(f"\n{target} percent identity")
        print(df.loc[df['seq1'] == 'Mycobacterium_intracellulare', 'Mycobacterium_chimaera'].values[0])
        data.append(df.loc[df['seq1'] == 'Mycobacterium_intracellulare', 'Mycobacterium_chimaera'].values[0])
    print(data)
    resDF.loc[1] = data
    print(resDF)
    resDF.to_csv("chimaera_vs_intracellulare.csv")


if __name__ == "__main__":
    main()