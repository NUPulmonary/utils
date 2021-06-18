__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os
from snakemake.shell import shell

input = snakemake.input
output = snakemake.output
cellranger_dir = snakemake.params.get("cellranger_dir", "")
transcriptome = snakemake.params.get("transcriptome", None)
chemistry = snakemake.params.get("chemistry", None)
sample = snakemake.wildcards.get("sample", None)

assert transcriptome is not None, "param transcriptome is required"
assert chemistry is not None, "param chemistry is required"
assert sample is not None, "sample wildcard is required"

input_paths = ",".join([os.path.realpath(i) for i in input])

load_cellranger = ""
if cellranger_dir == "":
    load_cellranger = "module load cellranger"
else:
    cellranger_dir += "/bin/"

# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    """
    module purge all
    {load_cellranger}

    outdir=`dirname "{output[0]}"`
    mkdir -p "$outdir"
    cd "$outdir"

    {cellranger_dir}cellranger count --id {sample} \
        --sample={sample} \
        --transcriptome={transcriptome} \
        --fastqs={input_paths} \
        --chemistry={chemistry} {log}
    """
)
