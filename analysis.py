import concatenate, snpDists, seqID, alignCheck, snpPlots

def main():
    # # make concatenated align
    print("concatenating alignment...")
    concat_targets = ['atpD', 'groL', 'rpoB', 'tuf', '16S', '23S']
    concatenate.main(concat_targets)

    # # run analysis on alignments
    print("calculating SNPs...")
    targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf', 'ITS', 'concat']
    snpDists.main(targets)

    print("calculating % identity...")
    seqID.main(targets)

    print("making plots...")
    alignCheck.main(targets)
    snpPlots.main(targets)

if __name__ == '__main__':
    main()