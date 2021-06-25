__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os

from snakemake.shell import shell


input = snakemake.input
output = snakemake.output
params = snakemake.params
names = params.get("names", {})
priorities = params.get("priorities", {})


assert len(output) == 1, "Expect single output as built cellbrowser folder"


# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

for i in input:
    sample = os.path.basename(i)
    conf_path = os.path.join(i, "cellbrowser.conf")
    conf = open(conf_path, "rt").read()
    name = names.get(sample, None)
    if name:
        conf = conf.replace(f"shortLabel='{sample}'", f"shortLabel='{name}'")
    priority = priorities.get(sample, None)
    if priority:
        conf += f"\npriority={priority}\n"
    open(conf_path, "wt").write(conf)

# TODO: move to subprocess or ensure input is parsed correctly
shell(
    """
    for i in {input}; do
        cbBuild -i $i/cellbrowser.conf -o {output} {log}
    done
    """
)
