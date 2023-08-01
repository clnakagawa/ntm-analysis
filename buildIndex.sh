#!/bin/bash

# arg1 for sp arg2 for db name (id)

mkdir "mmInds/$1"
# species name (e.g. Mycolicibacterium_fortuitum)
SP_NAME=$1
# shorter name (e.g. fortuitum)
SP_ID=$2

build_target_ind()
{
   minimap2 -x map-ont -d "mmInds/$SP_NAME/$SP_ID_$1.mmi" "mmRef/$SP_NAME/$1.txt"
}

for target in "16S" "23S" "atpD" "groL" "rpoB" "tuf"
do
   build_target_ind $target
done

