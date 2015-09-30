#!/bin/bash

#SBATCH --partition=serial_requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24000
#SBATCH --time=24:00:00
#SBATCH --output=signalp.%j.out
#SBATCH --error=signalp.%j.err
#SBATCH -J signalp

module load hpc/signalp-4.1

signalp -v ../compara2/all_prot.fa > all.signalp.out



