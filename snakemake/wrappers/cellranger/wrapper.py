__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os
import tempfile
import warnings
import pandas as pd
from snakemake.shell import shell

input = snakemake.input
output = snakemake.output
cellranger_dir = snakemake.params.get("cellranger_dir", "")
transcriptome = snakemake.params.get("transcriptome", None)
chemistry = snakemake.params.get("chemistry", None)
sample = snakemake.wildcards.get("sample", None)
skip_sample = snakemake.params.get("skip_sample", False)
sample_csv_path = snakemake.params.get("sample_csv_path", "")

input_fastq_type = snakemake.params.get("input_fastq_type", "gex")
gex_fastqs = snakemake.params.get("gex_fastqs", None)
antibodies = snakemake.params.get("antibodies", None)

#handle expected cells
#if no cell numbers are specified, default to 3000 to match cellranger defaults
#keeps backward compatibility with old versions as well
samples = pd.read_csv(sample_csv_path)
if "Expected" not in samples.columns:
    warnings.warn("No Expected column detected. Defaulting to 3000 cells/sample.")
    samples['Expected'] = 3000
expected_cells = samples[samples.Sample == sample, 'Expected']

# mode of count operation
# gex
# antibody
# gex+antibody
# see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis
mode = "gex"

if input_fastq_type == "antibody":
    assert antibodies is not None, "param antibodies is required for this mode of cellranger count"
    if gex_fastqs is None:
        mode = "antibody"
    else:
        mode = "gex+antibody"
elif input_fastq_type == "gex":
    pass
else:
    raise ValueError(f"Unrecognized input_fastq_type value `{input_fastq_type}', expected gex, antibody, gex+antibody")

assert transcriptome is not None, "param transcriptome is required"
assert chemistry is not None, "param chemistry is required"
assert sample is not None, "sample wildcard is required"

load_cellranger = ""
if cellranger_dir == "":
    load_cellranger = "module load cellranger"
else:
    cellranger_dir += "/bin/"

input_paths = [os.path.realpath(i) for i in input]
if gex_fastqs is not None:
    gex_fastqs = [os.path.realpath(i) for i in gex_fastqs]

sample_arg = ""
feature_ref = ""
lib_msg = ""

if mode == "gex":
    input_paths = ",".join(input_paths)
    input_arg = f"--fastqs={input_paths}"
    if not skip_sample:
        sample_arg = f"--sample={sample}"

if mode == "gex+antibody":
    antibodies = os.path.realpath(antibodies)
    n_fastq = len(input_paths) + len(gex_fastqs)
    libraries = pd.DataFrame({
        "fastqs": input_paths + gex_fastqs,
        "sample": sample,
        "library_type": [input_fastq_type] * len(input_paths) + ["gex"] * len(gex_fastqs),
    })
    libraries.library_type.replace({
        "gex": "Gene Expression",
        "antibody": "Antibody Capture",
    }, inplace=True)
    _, lib_path = tempfile.mkstemp()
    libraries.to_csv(lib_path, index=False)
    input_arg = f"--libraries={lib_path}"
    feature_ref = f"--feature-ref=\"{antibodies}\""
    lib_msg = f"echo \"Created libraries.csv at {lib_path}\""

# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    """
    module purge all
    {load_cellranger}

    echo "Running cellranger count in mode {mode}"
    {lib_msg}

    outdir=`dirname "{output[0]}"`
    mkdir -p "$outdir"
    cd "$outdir"

    {cellranger_dir}cellranger count --id {sample} \
        --transcriptome={transcriptome} \
        --chemistry={chemistry} \
        {sample_arg} \
        {input_arg} \
        {feature_ref} \
        {log}
    """
)
