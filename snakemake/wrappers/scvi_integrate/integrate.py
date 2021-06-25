import os

import pandas as pd
import scanpy as sc
import scvi

import sc_utils


def rename_genes(names):
    names = names.str.replace("^GRCh38_+", "", regex=True)
    names = names.str.replace("^SARS2_+", "SARS-CoV-2-", regex=True)
    names = names.str.replace("SARS-CoV-2-antisense", "Antisense")
    return names


def read_sample(
    sample_meta,
    min_genes=200,
    min_cells=5,
):
    sample = sample_meta.Sample
    sample_dir = sample_meta.__dir

    ds = sc.read_10x_h5(os.path.join(sample_dir, sample, "outs", "filtered_feature_bc_matrix.h5"))
    ds.var_names = rename_genes(ds.var_names)
    ds.var_names_make_unique(join=".")
    ds.obs_names = sample + "_" + ds.obs_names.str.replace("-\d$", "", regex=True)
    sc.pp.filter_cells(ds, min_genes=min_genes)
    sc.pp.filter_genes(ds, min_cells=min_cells)

    if "__filter_cells" in sample_meta and type(sample_meta.__filter_cells) == str:
        sample_cells = pd.read_csv(sample_meta.__filter_cells, index_col=0)
        sample_cells.index = sample_cells.index.str.replace("-\d+$", "", regex=True)
        ds = ds[ds.obs_names.isin(sample_cells.index), :]

    ds.layers["counts"] = ds.X
    for k, v in sample_meta.iteritems():
        if k.startswith("__"):
            continue
        ds.obs[k] = v
    return ds


def integrate_scvi(
    h5ad,
    meta,
    markers,
    samples_meta=None,
    batch_key=None,
    covariates=None,
    n_latent=50,
    max_epochs=400,
    model_path=None,
    resolution=0.5,
    scvi_kwargs=None
):
    if batch_key is None:
        raise ValueError("Please provide a batch_key")
    if covariates is None:
        covariates = []

    if samples_meta is None:
        raise ValueError("Please provide samples_meta with samples metadata for integration")
    samples_meta = pd.read_csv(samples_meta)

    if scvi_kwargs is None:
        scvi_kwargs = {}

    datasets = []
    for _, sample_meta in samples_meta.iterrows():
        datasets.append(read_sample(sample_meta))

    ds = datasets[0].concatenate(datasets[1:], join="outer")
    ds.var["mito"] = ds.var_names.str.startswith("MT-")
    ds.var["ribo"] = ds.var_names.str.match("^RP(L|S)")
    sc_utils.merge_gene_info(ds)
    sc.pp.calculate_qc_metrics(
        ds,
        qc_vars=["mito", "ribo"],
        percent_top=[10, 20],
        log1p=False,
        inplace=True
    )
    ds.obs[batch_key] = ds.obs[batch_key].astype(str)
    for k in covariates:
        ds.obs[k] = ds.obs[k].astype(str)

    # Setup and train scVI model
    scvi.data.setup_anndata(ds, layer="counts", batch_key=batch_key, categorical_covariate_keys=covariates)
    vae = scvi.model.SCVI(ds, n_latent=n_latent, **scvi_kwargs)
    vae.train(max_epochs=max_epochs)
    if model_path:
        vae.save(model_path, overwrite=True)

    # Get latent representation and process kNN, leiden & UMAP on them
    ds.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(ds, use_rep="X_scVI")
    sc.tl.leiden(ds, resolution=resolution)
    sc.tl.umap(ds)

    # Copy normalized counts to the dataset and put them to raw for export
    ds.X = vae.get_normalized_expression(
        library_size=10e4
    )
    sc.pp.log1p(ds)
    ds.raw = ds

    # Compute markers with scVI
    m = vae.differential_expression(groupby="leiden", mode="change", batch_correction=True)
    m["leiden"] = m.comparison.str.replace(" vs Rest", "")
    m = m.reset_index()
    m["gene"] = m["index"]

    m = m.loc[
        m["is_de_fdr_0.05"] & (m.lfc_mean > 0) & (m.non_zeros_proportion1 > 0.05),
        ["proba_not_de", "lfc_mean", "non_zeros_proportion1", "non_zeros_proportion2",
        "proba_not_de", "leiden", "gene"]
    ]

    m.columns = "p_val,avg_logFC,pct.1,pct.2,p_val_adj,cluster,gene".split(",")
    m = m.sort_values(["cluster", "avg_logFC"], ascending=[True, False])

    # Save object, metadata and markers
    ds.write_h5ad(h5ad)
    ds.obs.to_csv(meta)
    m.to_csv(markers)
