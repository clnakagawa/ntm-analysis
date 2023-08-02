#!/bin/bash

# arg1 for sp arg2 for db name (id)

mkdir "mmInds/$1"
# species name (e.g. Mycolicibacterium_fortuitum)
SP_NAME=$1
# shorter name (e.g. fortuitum)
SP_ID=$2

build_target_ind()
{
   minimap2 -x map-ont -d "mmInds/$SP_NAME/${SP_ID}_$1.mmi" "mmRef/$SP_NAME/$1.txt"
}


while IFS="," read -r rec1
do 
   echo $rec1
   build_target_ind $rec1
done < <(cut -d "," -f2 targets.csv | tail -n +2)
