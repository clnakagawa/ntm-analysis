#!/bin/bash

# arg1 for sp arg2 for db name (id)

mkdir "mmInds/$1"

./minimap2/minimap2 -x map-ont -d "mmInds/$1/$2_16S.mmi" "mmRef/$1/16S.txt"
./minimap2/minimap2 -x map-ont -d "mmInds/$1/$2_23S.mmi" "mmRef/$1/23S.txt"
./minimap2/minimap2 -x map-ont -d "mmInds/$1/$2_atpD.mmi" "mmRef/$1/atpD.txt"
./minimap2/minimap2 -x map-ont -d "mmInds/$1/$2_groL.mmi" "mmRef/$1/groL.txt"
./minimap2/minimap2 -x map-ont -d "mmInds/$1/$2_rpoB.mmi" "mmRef/$1/rpoB.txt"
./minimap2/minimap2 -x map-ont -d "mmInds/$1/$2_tuf.mmi" "mmRef/$1/tuf.txt"


