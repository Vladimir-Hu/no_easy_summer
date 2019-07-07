#!/usr/bin/bash
echo "==>Compiling programs to simulate gene expression."
g++ -O2 -fopenmp TAL_Gene.cpp -o tau_gene
g++ -O2 -fopenmp SSA_Gene.cpp -o ssa_gene
echo "==>Running simulation."
time ./tau_gene &> tau_gene.log
time ./ssa_gene &> ssa_gene.log
echo "Cleanning files."
mkdir -p ../data/Gene
mv raw_data_*.txt ../data/Gene
mv *_gene.log ../data/Gene