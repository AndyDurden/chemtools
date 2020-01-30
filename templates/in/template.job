#!/bin/bash
#SBATCH --job-name=tempname
#SBATCH --output=res_%j.txt
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16000
#SBATCH --gres gpu:1

module load intel
source /home/xsede/users/xs-adurden1/.bashrc



cd temppath

srun /cstor/xsede/users/xs-adurden1/terachem_xstream/terachem tempname.in > tempname.out

