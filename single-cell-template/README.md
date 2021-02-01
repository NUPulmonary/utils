
1. `module load python-anaconda3/2019.10`
2. `pipenv install`
3. `pipenv run snakemake --cluster /projects/b1038/tools/snakemake/cluster.py -j 50 --cluster-config ~/.slurm.json`
