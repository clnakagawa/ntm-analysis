import numpy as np
import pandas as pd

tbIDs = "tuberculosis|canettii|mungi|orygis"

def keyF(val):
    avals = ['Complete Genome', 'Chromosome']
    if val in avals:
        return avals.index(val)
    return 2


def getSp(x):
    x.replace("[|]","")
    strs = x.split(" ")
    splits = ['subsp.', 'variant']
    if len(strs) > 2 and strs[2] in splits:
        return "_".join([strs[0], strs[3]])
    return "_".join(strs[:2])

def main():
    pd.set_option('display.max_columns', 15)
    pd.set_option('display.max_rows', 10)
    assemblies = pd.read_table("assembly_summary_refseq.txt", low_memory=False, skiprows=1)
    data = assemblies[['# assembly_accession',
                       "taxid",
                       "species_taxid",
                       "organism_name",
                       "infraspecific_name",
                       "assembly_level",
                       "genome_rep",
                       "ftp_path"]]

    data = data.loc[data['genome_rep'] == "Full"]
    data = data.loc[data['organism_name'].str.startswith(("Mycobacterium","Mycobacteroides","Mycolicibacter","Mycolicibacillus","Mycolicibacterium"))]
    data['sp'] = data['organism_name'].apply(getSp)
    print(data.head())
    data = data.loc[data['assembly_level'] == 'Complete Genome']
    data['oKey'] = data['organism_name'].apply(lambda x: 1-('subsp' in x))
    print(data.head())
    data = data.sort_values(by=['oKey'])
    data = data.loc[data['organism_name'].str.contains("phage|Phage| sp. ") == False]
    data = data.loc[data['organism_name'].str.contains(tbIDs) == False]
    data = data.drop_duplicates(['sp'], keep="first")
    data.to_csv("summaryData.csv")
    print(f"# genomes: {len(data.index)}")

# pull sequences for testing
def test():
    pd.set_option('display.max_columns', 15)
    pd.set_option('display.max_rows', 10)
    assemblies = pd.read_table("assembly_summary_refseq.txt", low_memory=False, skiprows=1)
    data = assemblies[['# assembly_accession',
                       "taxid",
                       "species_taxid",
                       "organism_name",
                       "infraspecific_name",
                       "assembly_level",
                       "genome_rep",
                       "ftp_path"]]

    data = data.loc[data['genome_rep'] == "Full"]
    data = data.loc[data['organism_name'].str.startswith(("Mycobacterium","Mycobacteroides","Mycolicibacter","Mycolicibacillus","Mycolicibacterium"))]
    data['sp'] = data['organism_name'].apply(getSp)
    print(data.head())
    data = data.loc[data['assembly_level'] == 'Complete Genome']
    data['oKey'] = data['organism_name'].apply(lambda x: 1-('subsp' in x))
    print(data.head())
    data = data.sort_values(by=['oKey'])
    data = data.loc[data['organism_name'].str.contains("phage|Phage| sp. ") == False]
    data = data.loc[data['organism_name'].str.contains(tbIDs) == False]
    data.to_csv("testFiles/summaryData.csv")
    print(f"# genomes: {len(data.index)}")

if __name__ == "__main__":
    main()