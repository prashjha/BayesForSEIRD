#!/bin/bash

module load python3

# source conda
source /work/06743/prashjha/stampede2/Softwares/anaconda/etc/profile.d/conda.sh

# activate conda environment
conda activate fenics

python ft01_poisson.py
