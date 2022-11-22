#!/bin/bash

# This is a bash script used to print the SECOND k-vector (in cart coords.)
# in multiple CASINO out files. This is particularly useful in producing the twist
# vector list to be used in other scripts in this folder, e.g. merge_expval.py-script.

# In the following, I will assume that we have a set of folders containing CASINO
# simulations, ordered as i1,i1+1,...,i1+N and corresponding to different twists. 



if [ "$#" -ne 2 ]; then
    echo " "
    echo "Usage: "
    echo "./script <ind1> <ind2>, "
    echo "where"
    echo "<ind1>:      the index of the first CASINO folder,"
    echo "<ind2>:      the index of the last CASINO folder,"
    echo " "
    echo "CASINO out-files will be looked as ind/out "
    echo " "
    exit
fi

echo " "
echo "----- Produce a list of twists vectors -------"
echo " "

echo "Files to be parsed: "
for ((i=$1;i<=$2;i++))
do
    file=$i/out
    echo "- $file"
    grep "k(au)" $file >> f1
    head -2 f1 >> f2
    tail -1 f2 >> f3
    rm -f f1 f2
done

echo " "
echo "List of twist vectors: " 
awk '{print $2 " " $3 " " $4}' f3
rm -f f3
echo " "
echo "NOTE! The first k-vector should probably be checked, as gamma-twist is listed differently "
echo "      in the CASINO out-file."
echo " "
