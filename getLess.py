import numpy as np
import pandas as pd

def main():
    # targets = ['16S', '23S', 'atpD', 'groL', 'rpoB', 'tuf', 'ITS_allseqs']
    #
    # with open("primerSeqs/shortList.txt") as f:
    #     terms = f.read().split("\n")
    # print(terms)
    #
    # for target in targets:
    #     with open(f"lpsnFiles/cleanedSeqs/{target}.txt", 'r') as f:
    #         seqs = f.read().split("\n>")
    #     seqs = [seq for seq in seqs if any([term for term in terms if term in seq])]
    #     with open(f"primerSeqs/{target}.txt", 'w') as f:
    #         f.write(">"+"\n>".join(seqs))

    with open("primerSeqs/16S_align", 'r') as f:
        seqs16s = f.read().split("\n>")
    with open("primerSeqs/23S_align", 'r') as f:
        seqs23s = f.read().split("\n>")

    just16seqs = {seq.split("\n")[0].replace(">",''):''.join(seq.split("\n")[1:])for seq in seqs16s}
    just23seqs = {seq.split("\n")[0].replace(">",''):''.join(seq.split("\n")[1:]) for seq in seqs23s}
    names = [seq.split("\n")[0].replace(">",'') for seq in seqs16s]
    print(names)
    mod16seqs = [name+"\n"+just16seqs.get(name)+just23seqs.get(name)[:300] for name in names]
    mod23seqs = [name+"\n"+just23seqs.get(name)[300:] for name in names]
    print("\n>".join(mod16seqs))
    print("\n>".join(mod23seqs))
    with open("16Smod_align", 'w') as f:
        f.write(">"+"\n>".join(mod16seqs))
    with open("23mod_align", 'w') as f:
        f.write(">"+"\n>".join(mod23seqs))
if __name__ == "__main__":
    main()