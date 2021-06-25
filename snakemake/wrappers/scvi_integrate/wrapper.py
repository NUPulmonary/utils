__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os

import pandas as pd

import integrate


# For GPU on quest
# https://github.com/ray-project/ray/issues/10995#issuecomment-698177711
os.environ["SLURM_JOB_NAME"] = "bash"

input = snakemake.input
output = snakemake.output
params = snakemake.params
batch_key = params.get("batch_key", None)
covariates = params.get("covariates", [])
resolution = params.get("resolution", None)
model_path = params.get("model_path", None)
n_latent = params.get("n_latent", 50)
scvi_kwargs = params.get("scvi_kwargs", {})

assert len(input) == 1, "Expect single input of samples_meta csv file"
assert hasattr(output, "h5ad"), "Output should have `h5ad' field with path to h5ad file to save"
assert hasattr(output, "meta"), "Output should have `meta' field with path to metadata file to save"
assert hasattr(output, "markers"), "Output should have `markers' field with path to markers file to save"

assert resolution is not None, "Expect `resolution' parameter"
assert model_path is not None, "Expect `model_path' parameter, directory name where to save model"
assert model_path.endswith(".model"), "`model_path' should end with .model by convention"

integrate.integrate_scvi(
    **output,
    samples_meta=input[0],
    batch_key=batch_key,
    covariates=covariates,
    resolution=resolution,
    model_path=model_path,
    n_latent=n_latent,
    scvi_kwargs=scvi_kwargs,
)
