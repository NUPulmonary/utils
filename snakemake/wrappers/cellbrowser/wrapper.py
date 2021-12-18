__author__ = "Nikolay Markov"
__email__ = "nikolay.markov@northwestern.edu"
__license__ = "MIT"


import os

import pandas as pd
from snakemake.shell import shell


input = snakemake.input
output = snakemake.output
params = snakemake.params
sample = snakemake.wildcards.get("sample", None)

assert hasattr(input, "h5ad"), "Output should have `h5ad' field with path to h5ad file to save"
assert hasattr(input, "meta"), "Output should have `meta' field with path to metadata file to save"
assert hasattr(input, "markers"), "Output should have `markers' field with path to markers file to save"
assert len(output) == 1, "Expect single output with tar.gz file"

assert sample is not None, "sample wildcard is required"

cluster_field = params.get("cluster_field", "leiden")

# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    """
    h5=`realpath {input.h5ad}`
    out_file=`realpath {output}`
    base_dir=`dirname $out_file`

    meta=`realpath "{input.meta}"`
    meta_file=`basename $meta`
    markers=`realpath "{input.markers}"`
    markers_file=`basename $markers`

    mkdir -p $base_dir
    fname=`basename {output}`
    out=${{fname%.tar.gz}}
    out_dir="$base_dir/$out"
    rm -rf "$out_dir"

    cbImportScanpy -i $h5 -o $out_dir --skipMarkers {log}

    echo "Sedding " "$out_dir/cellbrowser.conf"
    sed -i 's/louvain/leiden/g' "$out_dir/cellbrowser.conf"
    sed -i "s/defColorField='leiden'/defColorField='{cluster_field}'/" "$out_dir/cellbrowser.conf"
    sed -i "s/labelField='leiden'/labelField='{cluster_field}'/" "$out_dir/cellbrowser.conf"
    sed -i 's/#radius=2/radius=2/' "$out_dir/cellbrowser.conf"
    sed -i "s/meta.tsv/$meta_file/" "$out_dir/cellbrowser.conf"
    sed -i -z 's/ \+{{\\n \+"file": "scVI_coords.tsv"[^}}]\+}},\\n//' "$out_dir/cellbrowser.conf"
    sed -i -z 's/ \+{{\\n \+"file": "_scvi[^}}]\+}},\?\\n//' "$out_dir/cellbrowser.conf"

    echo -e "\ndisplay_single_meta=True" >> "$out_dir/cellbrowser.conf"
    echo -e "\nmarkers = [{{\\"file\\": \\"$markers_file\\", \\"shortLabel\\":\\"Cluster Markers\\"}}]\n" \
        >> "$out_dir/cellbrowser.conf"

    cp "$meta" "$markers" $out_dir
    """
)
