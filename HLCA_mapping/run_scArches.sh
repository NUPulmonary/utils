#!/bin/bash

#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --mem=150G
#SBATCH --ntasks-per-node=12
#SBATCH --gres=gpu:a100:1
#SBATCH --mail-user=samuel.fenske@northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name='new_hlca'
#SBATCH --output=new_hlca.out
#SBATCH --error=new_hlca.err

script_dir=$(pwd)
code_dir=$'/projects/b1038/Pulmonary/sfenske/projects/scArches/code'

#parameters to edit
input_h5=$'/projects/b1038/Pulmonary/sfenske/sequencing/data/2022-44/cellranger/SC503/outs/filtered_feature_bc_matrix.h5'
dir_out=$'/projects/b1038/Pulmonary/sfenske/projects/scArches/output'
sample='SC503'
overwrite_models='yes' #overwrite model for sample if model already exists ('no' or 'yes')

#activate conda environment
cd $code_dir
eval "$(conda shell.bash hook)"
conda activate HLCA_mapping_env

#run python script
cd script_dir
python3 scArches.py --path $input_h5 --sample $sample --overwrite_models $overwrite_models  --dir_out $dir_out > scArches_log.txt

