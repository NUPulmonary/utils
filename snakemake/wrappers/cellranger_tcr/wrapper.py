__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os
import tempfile

import pandas as pd
from snakemake.shell import shell

input = snakemake.input
output = snakemake.output
cellranger_dir = snakemake.params.get("cellranger_dir", "")
transcriptome = snakemake.params.get("transcriptome", None)
sample = snakemake.wildcards.get("sample", None)

assert transcriptome is not None, "param transcriptome is required"
assert sample is not None, "sample wildcard is required"

load_cellranger = ""
if cellranger_dir == "":
    load_cellranger = "module load cellranger"
else:
    cellranger_dir += "/bin/"

input_paths = [os.path.realpath(i) for i in input]

feature_ref = ""

input_paths = ",".join(input_paths)
input_arg = f"--fastqs={input_paths}"

# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    """
    module purge all
    {load_cellranger}

    outdir=`dirname "{output[0]}"`
    mkdir -p "$outdir"
    cd "$outdir"

    {cellranger_dir}cellranger vdj --id {sample} \
        --reference={transcriptome} \
        --sample={sample} \
        {input_arg} \
        {log}
    """
)
