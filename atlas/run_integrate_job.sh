#!/bin/bash
#SBATCH -J scvi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tesla:1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=164GB

ml any/python/3.9.9
source /home/scvi_venv/bin/activate
cd /home/immune_ageing/HPC/

python integrate.py -i /home/immune_ageing_in/prepared.h5ad -o /home/immune_ageing_out
