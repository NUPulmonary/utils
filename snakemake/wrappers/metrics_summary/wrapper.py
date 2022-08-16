__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os
import pandas as pd

input = snakemake.input
output = snakemake.output

result = []
for sample in input:
    m = pd.read_csv(f"{sample}/outs/metrics_summary.csv")
    m.insert(0, "Sample", os.path.basename(sample))
    m.insert(2,'Fraction Reads in Cells',m.pop('Fraction Reads in Cells'))
    result.append(m)
result = sorted(result, key=lambda x: x.shape[1], reverse=True)
result = pd.concat(result, sort=False)
for column in result.columns:
    if result[column].dtype not in ("int64", "float64"):
        match = (~result[column].isna()) & result[column].str.match(r"^(\d+,)*\d+$")
        new_column = result[column].copy()
        new_column[match] = new_column.loc[match].str.replace(",", "")
        result[column] = new_column
result.sort_values("Sample").to_csv(output[0], index=False)
