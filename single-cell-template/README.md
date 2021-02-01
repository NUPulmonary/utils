
1. Edit & copy `.slurm.json` to your home folder
2. `module load python-anaconda3/2019.10`
3. `pipenv install`: I recommend setting `PIPENV_VENV_IN_PROJECT=1` to have venv in the project folder
4. Edit the `samples.csv`
5. Run snakemake from venv `pipenv run snakemake --cluster /projects/b1038/tools/snakemake/cluster.py -j 50 --cluster-config ~/.slurm.json`
