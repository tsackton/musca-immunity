#!/bin/bash

#SBATCH -p general
#SBATCH --mem 8000
#SBATCH --time 36:00:00
#SBATCH -J make_ultra
#SBATCH -n 16
#SBATCH -o raxml_%j.out
#SBATCH -e raxml_%j.err

trimal -in ultrametric.fa -out ultrametric_raxml.phy -phylip -nogaps
raxmlHPC-PTHREADS-AVX -n ULTRANG -f e -s ultrametric_raxml.phy -m PROTGAMMAAUTO -t dipt.stree -T 16 -p 12345 -N 10


