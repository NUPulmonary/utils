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

assert transcriptome is not None, "param transcriptome is required"
assert chemistry is not None, "param chemistry is required"

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

    {cellranger_dir}cellranger count --id {wildcards.sample} \
        --sample={wildcards.sample} \
        --transcriptome={params.transcriptome} \
        --fastqs={params.input_paths} \
        --chemistry={params.chemistry} {log}
    """
)
