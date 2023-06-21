import numpy as np
import pandas as pd
import getGenes, extractGene, os

def getSp(x):
    x.replace("[","")
    x.replace("]","")
    strs = x.split(" ")
    splits = ['subsp.', 'variant']
    if len(strs) > 2 and strs[2] in splits:
        return "_".join([strs[0], strs[3]])
    return "_".join(strs[:2])

def keyF(val):
    avals = ['Complete Genome', 'Chromosome']
    if val in avals:
        return avals.index(val)
    return 2

def spName(row):
    if pd.isnull(row['subsp_epithet']):
        return f"{row['genus_name']} {row['sp_epithet']}"
    return f"{row['genus_name']} {row['subsp_epithet']}"

def main():
    # find ntms in lsbn csv
    # get strain ids for type strains
    ntmGenus = ['Mycobacterium', 'Mycobacteroides', 'Mycolicibacter', 'Mycolicibacterium', 'Mycolicibacillus']
    tbNames = [' tuberculosis', ' africanum', ' mungi', ' orygis', ' microti', ' pinnipedii', ' bovis']
    pd.set_option('display.max_columns', 20)
    fname = "lpsn_gss_2023-06-20.csv"
    lpsn = pd.read_csv(fname)
    lpsn = lpsn.loc[lpsn['genus_name'].str.contains("|".join(ntmGenus))]
    data = [[spName(row),
             '|'.join(row['nomenclatural_type'].split("; ")),
             pd.isnull(row['subsp_epithet'])] for ind, row in lpsn.iterrows()]
    typeStrains = pd.DataFrame(data, columns=['spName', 'IDs', 'noSubsp'])
    typeStrains = typeStrains.sort_values(by='noSubsp')
    typeStrains = typeStrains.drop_duplicates('spName')
    typeStrains = typeStrains.loc[~typeStrains['spName'].str.contains('|'.join(tbNames))]
    typeStrains.to_csv("lbsnNTM.csv")

    # find refseq entries for typestrains if exist
    refseqs = pd.read_table("assembly_summary_refseq.txt", low_memory=False, skiprows=1)
    refseqs = refseqs.loc[refseqs['organism_name'].str.contains('|'.join(ntmGenus))]
    idx = [refseqs.loc[(refseqs['infraspecific_name'].str.contains(row['IDs'], na=False))
           & (refseqs['organism_name'].str.contains(row['spName'].split(' ')[-1]))].values for ind, row in typeStrains.iterrows()]
    idx = [item for sublist in idx for item in sublist]
    refTypes = pd.DataFrame(idx, columns=refseqs.columns)
    refTypes = refTypes.loc[~refTypes['organism_name'].str.contains('|'.join(tbNames))]
    refTypes['assembly_level'] = pd.Categorical(refTypes['assembly_level'], ['Complete Genome', 'Chromosome', 'Contig', 'Scaffold'])
    refTypes = refTypes.sort_values('assembly_level')
    refTypes['sp'] = refTypes['organism_name'].apply(getSp)
    refTypes = refTypes.drop_duplicates('taxid')
    refTypes.to_csv("refseqTypes.csv")

    if not os.path.exists("lbsnFiles"):
        os.makedirs("lbsnFiles")

    # pull files
    getGenes.main(file="refseqTypes.csv", output="lbsnFiles", type="cds")
    getGenes.main(file="refseqTypes.csv", output="lbsnFiles", type="rna")
    getGenes.main(file="refseqTypes.csv", output="lbsnFiles")

    # extract targets
    targets = [('chaperonin GroEL', 'groL'),
               ('ATP synthase subunit beta', 'atpD'),
               ('RNA polymerase subunit beta', 'rpoB'),
               ('elongation factor Tu', 'tuf')]
    extractGene.main(targets, file="refseqTypes.csv", dir="lbsnFiles")





if __name__ == "__main__":
    main()