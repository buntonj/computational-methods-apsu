#! /bin/bash
#SBATCH -j TITLE
#SBATCH -o ~/LOCATION/OUTPUT.txt
#SBATCH -n 1
#SBATCH -p aux
#SBATCH --time=DD-HH:MM:SS

cd ~/LOCATION/CODE.whatever
