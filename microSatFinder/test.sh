#!/bin/sh
#SBATCH --job-name=test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=imir@ucsc.edu
#SBATCH --ntasks=1
#SBATCH --mem=600mb
#SBATCH --output=serial_test_%j.out

python3 findMicroSat.py -infile cenX_chm13.alphaSatDb.guppy_3.1.5.FR.rc.fa.gz
 
