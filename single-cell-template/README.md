
1. Edit & copy `.slurm.json` to your home folder
2. `module load python-anaconda3/2019.10`
3. Ensure `pipenv` is installed: `pip install --user pipenv`
4. Copy files to your project folder. Next items should be executed in the project folder
5. `pipenv install`: I recommend using `PIPENV_VENV_IN_PROJECT=1 pipenv install` to have venv in the project folder
6. Edit the `samples.csv`
7. Run snakemake from venv `pipenv run snakemake --cluster /projects/b1038/tools/snakemake/cluster.py -j 50 --cluster-config ~/.slurm.json`
