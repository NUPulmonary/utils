__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import pandas as pd

input = snakemake.input
output = snakemake.output
run_id = snakemake.wildcards.get("run_id", None)
samples = pd.read_csv(input[0])

if run_id is None:
    raise ValueError("Expect run_id wildcard in output path")

columns = ["Lane", "Sample", "Index"]
full_columns = pd.Series(["RunID"] + columns)
missing_columns = full_columns[~full_columns.isin(samples.columns)].values

if missing_columns.size > 0:
    raise ValueError(f"Could not find columns {', '.join(missing_columns)} in {input}")

samples = samples.loc[samples.RunID == run_id, columns]
samples.columns = columns
samples.to_csv(output[0], index=False)
