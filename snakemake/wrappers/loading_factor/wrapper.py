__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os
import pandas as pd

input = snakemake.input
output = snakemake.output

data = pd.read_csv(input[0])
data = data.iloc[:, :4]
depth = data["Mean Reads per Cell"]
data["Loading factor"] = depth.max() / depth
data.to_csv(output[0], index=False)
