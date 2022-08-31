__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os
from snakemake.shell import shell

input = snakemake.input
output = snakemake.output
run_id = snakemake.wildcards.get("run_id", None)
cellranger_dir = snakemake.params.get("cellranger_dir", "")
rc_i2_override = snakemake.params.get("r2_i2_override", None)
keep_output = snakemake.params.get("keep_cellranger_output", False)

assert len(input) == 2, "expecting 2 inputs: bcl directory and sample sheet"
assert run_id is not None, "run_id is a required wildcard"

flowcell = run_id.split("_")[-1][1:]

load_cellranger = ""
if cellranger_dir == "":
    load_cellranger = "module load cellranger"
else:
    cellranger_dir += "/bin/"

rm_output = ""
if not keep_output:
    rm_output = f"rm -rf {flowcell}\nrm -f __{flowcell}.mro"

# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# check if this is MiniSeq and we need to set rc-i2-override
# because mkfastq does not set it
if rc_i2_override is None:
    rc_i2_override=""
    for l in open(os.path.join(input[0], "RunParameters.xml")):
        if "ApplicationName" in l:
            if "MiniSeq" in l:
                rc_i2_override="--rc-i2-override=true"
            break
else:
    rc_i2_override=f"--rc-i2-override={rc_i2_override}"


# TODO: account for lanes
# TODO: merge lanes

shell(
    """
    module purge all
    module load bcl2fastq/2.20
    {load_cellranger}

    if [[ -d {flowcell} ]]; then
        echo "{flowcell} directory exists, removing it"
        rm -rf {flowcell}
        rm -f __{flowcell}.mro
    fi

    {cellranger_dir}cellranger mkfastq --run="{input[0]}" \
        --csv="{input[1]}" \
        {rc_i2_override} \
        --output-dir="{output[0]}" {log}

    {rm_output}
    """
)
