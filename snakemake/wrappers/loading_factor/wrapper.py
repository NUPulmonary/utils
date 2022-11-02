__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os
import pandas as pd

input = snakemake.input
output = snakemake.output
data = pd.read_csv(input[0])

#add handling for TCR and GEX metrics_summary.csv formats
if "Mean Reads per Cell" in data.columns:
    field = snakemake.params.get("field", "Mean Reads per Cell")
elif "Mean Read Pairs per Cell" in data.columns:
    field = snakemake.params.get("field", "Mean Read Pairs per Cell")

lf_data = data.iloc[:, :2]
lf_data[field] = data[field]
depth = data[field]
lf_data["Loading factor"] = depth.max() / depth
lf_data.to_csv(output[0], index=False)
