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
sample_types = snakemake.params.get("sample_types", None)
skip_sample = snakemake.params.get("skip_sample", False)
sample_csv_path = snakemake.params.get("sample_csv_path", None)
expected_cells = snakemake.params.get("n_cells", None)
include_introns = snakemake.params.get("include_introns", None)

input_fastq_type = snakemake.params.get("input_fastq_type", "gex")
gex_fastqs = snakemake.params.get("gex_fastqs", None)
antibodies = snakemake.params.get("antibodies", None)

#handle expected cells
if expected_cells is None and sample_csv_path is not None:
    samples = pd.read_csv(sample_csv_path)
    if 'Expected' in samples:
        expected_cells = samples[samples.Sample == sample, 'Expected'].values[0]

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
elif input_fastq_type == "gex+antibody":
    assert antibodies is not None, "param antibodies is required for this mode of cellranger count"
    assert sample_types is not None, "param sample_types is required for this model of cellranger count"
    mode = "gex+antibody"
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
    # We can invoke this with fastq files in different directories and same name
    # or with different names and same directory
    if gex_fastqs is None:  # same directory option
        libraries = pd.DataFrame({
            "fastqs": input_paths[0],  # ugly hack, fix later
            "sample": list(sample_types.values()),
            "library_type": list(sample_types.keys())
        })
    else:  # same name option
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
    antibodies = os.path.realpath(antibodies)
    print(libraries);
    input_arg = f"--libraries={lib_path}"
    feature_ref = f"--feature-ref=\"{antibodies}\""
    lib_msg = f"echo \"Created libraries.csv at {lib_path}\""

expected_cells_arg = ""
if expected_cells is not None:
    expected_cells_arg = f"--expect-cells {expected_cells}"

include_introns_arg = ""
if include_introns is not None:
    include_introns_arg = "--include-introns false"
    if include_introns:
        include_introns_arg = "--include-introns true"

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
        {include_introns_arg} \
        {sample_arg} \
        {input_arg} \
        {feature_ref} \
        {expected_cells_arg} \
        {log}
    """
)
