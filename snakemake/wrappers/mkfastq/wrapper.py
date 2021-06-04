__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os
from snakemake.shell import shell

input = snakemake.input
output = snakemake.output
run_id = snakemake.wildcards.get("run_id", None)
cellranger_dir = snakemake.params.get("cellranger_dir", "")

assert len(input) == 2, "expecting 2 inputs: bcl directory and sample sheet"
assert run_id is not None, "run_id is a required wildcard"

flowcell = run_id.split("_")[-1][1:]

load_cellranger = ""
if cellranger_dir == "":
    load_cellranger = "module load cellranger"
else:
    cellranger_dir += "/bin/"

# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# TODO: account for lanes
# TODO: check rc-i2-override

shell(
    """
    module purge all
    module load bcl2fastq/2.19.1
    {load_cellranger}

    {cellranger_dir}cellranger mkfastq --run="{input[0]}" \
        --csv="{input[1]}" \
        --rc-i2-override=true \
        --output-dir="{output[0]}" {log}

    rm -rf {flowcell}
    rm __{flowcell}.mro
    """
)
