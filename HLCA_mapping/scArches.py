import os
import scarches as sca
import scanpy as sc
import anndata
import sys
import gdown
import torch
import matplotlib.pyplot as plt
from scarches.dataset.trvae.data_handling import remove_sparsity
import numpy as np
import scipy
import pandas as pd
import scvi
import urllib.request 
import gzip
import shutil
import warnings
import argparse
import math

import sys
os.chdir('/projects/b1038/Pulmonary/sfenske/projects/scArches/code')
sys.path.append('../mapping_data_to_the_HLCA/scripts')
import scarches_label_transfer
import data_import_and_cleaning

parser=argparse.ArgumentParser()
parser.add_argument('--path',type=str,required=True)
parser.add_argument('--sample',type=str,required=True)
parser.add_argument('--overwrite_models',type=str,required=True)
parser.add_argument('--dir_out',type=str,required=True)

args=parser.parse_args()


path_gene_order = "../mapping_data_to_the_HLCA/supporting_files/HLCA_scarches_gene_order.csv"
path_embedding = "../mapping_data_to_the_HLCA/data/HLCA_emb_and_metadata.h5ad"
path_HLCA_celltype_info = "../mapping_data_to_the_HLCA/supporting_files/HLCA_celltypes_ordered.csv"
dir_ref_model = "../mapping_data_to_the_HLCA/data/HLCA_reference_model"

testdata = args.path
sample=args.sample
if args.overwrite_models=='yes':
    overwrite_models=True
else:
    overwrite_models=False

dir_out=args.dir_out
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

# gene order for scArches model
reference_gene_order = pd.read_csv(path_gene_order)
# reference embedding, including cell/sample/subject metadata:
reference_embedding = sc.read_h5ad(path_embedding)

def import_testdata(testdata):
    # read in adata
    
    if testdata.endswith('.h5'):
        query_data_full = sc.read_10x_h5(testdata)
        query_data_full.var_names_make_unique()
    
        # clean up .var.index (gene names)
        query_data_full.var['gene_names'] = query_data_full.var.index
        query_data_full.var.index = [idx.split("___")[-1] for idx in query_data_full.var.gene_ids]
        # clean up cell barcodes:
        query_data_full.obs.index = query_data_full.obs.index.str.rstrip("-1")
    
    elif testdata.endswith('.h5ad'):
        query_data_full = sc.read_h5ad(testdata)
    
    else:
        print('Expects an .h5 or .h5ad file as input')
        sys.exit()
    
    
    
#     # read in metadata (to select only cells of interest and remove empty drops)
#     query_data_meta = pd.read_csv(os.path.join(testfolder,"testmeta.csv.gz"),index_col=0)
#     # subset to cells from our sample
#     query_data_meta = query_data_meta.loc[query_data_meta.donor == "D12_4",:].copy()
#     # clean up barcodes:
#     query_data_meta.index = [idx.split("-")[-1] for idx in query_data_meta.index]
#     # subset adata to cells in metadata:
#     query_data_full = query_data_full[query_data_meta.index,:].copy()
    
    
    
    # add dataset information:
    # query_data_full.obs['dataset'] = "test_dataset_delorey_regev"
    
    
    clean_genes=[]
    for idx in query_data_full.var.index:
        if 'GRCh38' in idx:
            clean_genes.append(idx.replace('GRCh38_',''))
        else:
            clean_genes.append(idx)
    query_data_full.var.index=clean_genes
    
    query_data_full.var.gene_ids=query_data_full.var.gene_names

    query_data_full.var.gene_ids=query_data_full.var.gene_names

    return query_data_full
    
def subset_and_pad_adata_object(adata, ref_genes_df, min_n_genes_included=1500):
    adata.var=adata.raw.var
    # test if adata.var.index has gene names or ensembl names:
    n_ids_detected = sum(adata.var.index.isin(ref_genes_df.gene_id))
    n_symbols_detected = sum(adata.var.index.isin(ref_genes_df.gene_symbol))
    if max(n_ids_detected, n_symbols_detected) < min_n_genes_included:
        # change column names to lower case
        adata.var.columns = adata.var.columns.str.lower()
        # check if gene names are in another column:
        if "gene_symbols" in adata.var.columns:
            adata.var.index = adata.var.gene_symbol
            n_symbols_detected = sum(adata.var.index.isin(ref_genes_df.gene_symbol))
        elif "gene_ids" in adata.var.columns:
            adata.var.index = adata.var.gene_ids
            n_ids_detected = sum(adata.var.index.isin(ref_genes_df.gene_id))
        # check if anything changed:
        if max(n_ids_detected, n_symbols_detected) < min_n_genes_included:    
            raise ValueError(f"We could detect only {max(n_ids_detected, n_symbols_detected)} genes of the 2000 that we need for the mapping! The minimum overlap is {min_n_genes_included}. Contact the HLCA team for questions. Exiting")
    else:
        if n_ids_detected >= n_symbols_detected:
            gene_type = "gene_id"
            print("Gene names detected: ensembl gene ids.")
            n_genes = n_ids_detected  
        else:
            gene_type = "gene_symbol"
            n_genes = n_symbols_detected
            print("Gene names detected: ensembl gene symbols.")
    genes = adata.var.index[adata.var.index.isin(ref_genes_df[gene_type])].tolist()
    # if not all genes are included, pad:
    if n_genes > 2000:
        raise ValueError("Your gene names appear not to be unique, something must be wrong. Exiting.")
    print(f"{n_genes} genes detected out of 2000 used for mapping.")
    # Subset adata object
    adata_sub = adata[:,genes].copy()
    # Pad object with 0 genes if needed
    if n_genes < 2000:
        diff = 2000 - n_genes
        print(f'Not all genes were recovered, filling in zeros for {diff} missing genes...')
        # Genes to pad with
        genes_to_add = set(ref_genes_df[gene_type].values).difference(set(adata_sub.var_names))
        df_padding = pd.DataFrame(data=np.zeros((adata_sub.shape[0],len(genes_to_add))), index=adata_sub.obs_names, columns=genes_to_add)
        adata_padding = sc.AnnData(df_padding)
        # Concatenate object
        adata_sub = anndata.concat([adata_sub, adata_padding], axis=1, join='outer', index_unique=None, merge='unique')
        # and order:
        adata_sub = adata_sub[:,ref_genes_df[gene_type]].copy()
    return adata_sub
    
query_data_full = import_testdata(testdata) 
query_data_full.obs['dataset']=sample
query_data = data_import_and_cleaning.subset_and_pad_adata_object(query_data_full, reference_gene_order)
query_data.var=query_data.var.reindex(reference_gene_order.gene_id)
if (query_data.var.index == reference_gene_order.gene_symbol).all() or (query_data.var.index == reference_gene_order.gene_id).all():
    print("Gene order is correct.")
else:
    print(
        "WARNING: your gene order does not match the order of the HLCA reference. Fix this before continuing!"
    )
    
query_data.raw = query_data
raw = query_data.raw.to_adata()
raw.X = query_data.X
query_data.raw = raw
query_data.obs["scanvi_label"] = "unlabeled"
query_data.obs_names_make_unique()

batch_variable = "dataset"  # the column name under which you stored your batch variable
query_batches = sorted(query_data.obs[batch_variable].unique())

surgery_epochs = 500
early_stopping_kwargs_surgery = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "on": "full_dataset",
    "patience": 10,
    "threshold": 0.001,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

for batch in query_batches: # this loop is only necessary if you have multiple batches, but will also work for a single batch.
    print("Batch:", batch)
    surgery_path = os.path.join(dir_out,f"models/{batch}")
    if os.path.exists(surgery_path) and not overwrite_models:
        continue
    else:
        query_subadata = query_data[query_data.obs[batch_variable] == batch,:].copy()
        print("Shape:", query_subadata.shape)
        # load model and set relevant variables:
        model = sca.models.SCANVI.load_query_data(
            query_subadata,
            dir_ref_model,
            freeze_dropout = True,
            use_cuda=True
        )
        model._unlabeled_indices = np.arange(query_subadata.n_obs)
        model._labeled_indices = []
        # now train surgery model using reference model and target adata
        model.train(
            n_epochs_semisupervised=surgery_epochs,
            train_base_model=False,
            semisupervised_trainer_kwargs=dict(
                metrics_to_monitor=["accuracy", "elbo"], 
                weight_decay=0,
                early_stopping_kwargs=early_stopping_kwargs_surgery
            ),
            frequency=1
        )
        if not os.path.exists(surgery_path):
            os.makedirs(surgery_path)
        model.save(surgery_path, overwrite=True)
    
emb_df = pd.DataFrame(index=query_data.obs.index,columns=range(0,reference_embedding.shape[1]))

for batch in query_batches: # from small to large datasets
    print(f"Working on {batch}...")
    query_subadata = query_data[query_data.obs[batch_variable] == batch,:].copy()
    surgery_path = os.path.join(dir_out, f"models/{batch}")
    model = sca.models.SCANVI.load(surgery_path, query_subadata)
    query_subadata_latent = sc.AnnData(model.get_latent_representation(adata=query_subadata))
    # copy over .obs
    query_subadata_latent.obs = query_data.obs.loc[query_subadata.obs.index,:]
    query_subadata_latent.write(os.path.join(surgery_path, "emb.h5ad"))
    emb_df.loc[query_subadata.obs.index,:] = query_subadata_latent.X
    
query_embedding = sc.AnnData(X=emb_df.values, obs=query_data.obs)
query_embedding.obs['HLCA_or_query'] = "query"
reference_embedding.obs['HLCA_or_query'] = "HLCA"
#try altering reference matrix to latent space altogether...this is how old HCLA was setup
# reference_embedding=sc.AnnData(X=reference_embedding.obsm['X_scanvi_emb'], obs=reference_embedding.obs)

# combined_emb = query_embedding.concatenate(sc.AnnData(X=reference_embedding.obsm['X_scanvi_emb'], obs=reference_embedding.obs), index_unique=None) # index_unique="_", batch_key="dataset") # alternative
combined_emb = reference_embedding.concatenate(query_embedding, index_unique=None) # index_unique="_", batch_key="dataset") # alternative
combined_emb.X=np.nan_to_num(combined_emb.X)

sc.pp.neighbors(combined_emb, n_neighbors=30)
sc.tl.umap(combined_emb)

cts_ordered = pd.read_csv(path_HLCA_celltype_info,index_col=0)    

#clean up nan values before running weighted knn transfer
# ref_adata_obs = reference_embedding.obs.join(cts_ordered, on='ann_finest_level')
# for col in ref_adata_obs.columns[ref_adata_obs.columns.str.startswith('Level')]:
    # ref_adata_obs[col] = ['nan' if isinstance(x,float) and math.isnan(x) else x for x in ref_adata_obs[col]]
    
    
# run k-neighbors transformer
k_neighbors_transformer = scarches_label_transfer.weighted_knn_trainer(
    train_adata=reference_embedding,
    train_adata_emb="X", # location of our joint embedding
    n_neighbors=50,
    )    
# perform label transfer
labels, uncert = scarches_label_transfer.weighted_knn_transfer(
    query_adata=query_embedding,
    query_adata_emb="X", # location of our joint embedding
    label_keys="Level",
    knn_model=k_neighbors_transformer,
    ref_adata_obs = reference_embedding.obs.join(cts_ordered, on='ann_finest_level')
    )
    
uncertainty_threshold = 0.2
labels.rename(columns={f"Level_{lev}":f"Level_{lev}_transfered_label_unfiltered" for lev in range(1,6)},inplace=True)
uncert.rename(columns={f"Level_{lev}":f"Level_{lev}_transfer_uncert" for lev in range(1,6)},inplace=True)

combined_emb.obs = combined_emb.obs.join(labels)
combined_emb.obs = combined_emb.obs.join(uncert)
# convert to arrays instead of categoricals, and set "nan" to NaN 
for col in combined_emb.obs.columns:
    if col.endswith("_transfer_uncert"):
        combined_emb.obs[col] = list(np.array(combined_emb.obs[col]))
    elif col.endswith("_transfered_label_unfiltered"):
        filtered_colname = col.replace("_unfiltered","")
        matching_uncert_col = col.replace("transfered_label_unfiltered","transfer_uncert")
        
        # also create a filtered version, setting cells with too high 
        # uncertainty levels to "Unknown"
        combined_emb.obs[filtered_colname] = combined_emb.obs[col]
        combined_emb.obs[filtered_colname].loc[combined_emb.obs[matching_uncert_col]>uncertainty_threshold] = "Unknown"
        # convert to categorical:
        combined_emb.obs[col] = pd.Categorical(combined_emb.obs[col])
        combined_emb.obs[filtered_colname] = pd.Categorical(combined_emb.obs[filtered_colname])
        # then replace "nan" with NaN (that makes colors better in umap)
        combined_emb.obs[col].replace("nan",np.nan,inplace=True)
        combined_emb.obs[filtered_colname].replace("nan",np.nan,inplace=True)

print(f"Percentage of unknown per level, with uncertainty_threshold={uncertainty_threshold}:")
for level in range(1,6):
    print(f"Level {level}: {np.round(sum(combined_emb.obs[f'Level_{level}_transfered_label']=='Unknown')/query_data.n_obs*100,2)}%")

if not os.path.exists(f"{dir_out}/objects"):
    os.makedirs(f"{dir_out}/objects")
filtered_emb=combined_emb[combined_emb.obs.HLCA_or_query=='query']

# address following error
# Above error raised while writing key 'predicted_doublets' of <class 'h5py._hl.group.Group'> to /
for col in filtered_emb.obs.columns:
    if filtered_emb.obs[col].dtype=='category' or filtered_emb.obs[col].dtype=='object':
        filtered_emb.obs[col]=filtered_emb.obs[col].astype('str')
        
filtered_emb.write_h5ad(f"{dir_out}/objects/{sample}_emb.h5ad")


