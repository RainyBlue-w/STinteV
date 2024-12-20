import spatialdm as sdm
import anndata
import squidpy as sq
import pandas as pd

def calc_spatialDM(
    adata: anndata.AnnData, 
    species: str, l=20, cutoff=0.1, 
    n_nearest_neighbors=6, 
    single_cell=True, 
    min_cell=3,
    nproc = 1
):
    sdm.weight_matrix(
        adata, l=l, embedding='spatial', cutoff=cutoff, 
        n_nearest_neighbors=n_nearest_neighbors, single_cell=single_cell
    )
    sdm.extract_lr(adata, species=species, min_cell=min_cell)
    sdm.spatialdm_global(adata, n_perm=1000, method='z-score', nproc=nproc)
    sdm.sig_pairs(adata, method='z-score', fdr=True, threshold=0.05)
    sdm.spatialdm_local(adata, n_perm=1000, method='z-score', specified_ind=None, nproc=nproc)
    
    stats = {}
    for pair in adata.uns['global_res'][adata.uns['global_res'].fdr<0.05].sort_values(by='fdr').index :
        adata.obs['sq_autocorr'] = adata.uns['selected_spots'].loc[pair,:]
        stats[pair] = sq.gr.spatial_autocorr(adata, attr='obs', genes='sq_autocorr', connectivity_key='nearest_neighbors', mode='geary', copy=True)
        stats[pair]['%all'] = adata.uns['selected_spots'].loc[pair,:].mean() * 100
        stats[pair]['z'] = adata.uns['global_res'].loc[pair,'z']
        stats[pair]['fdr'] = adata.uns['global_res'].loc[pair,'fdr']
        stats = pd.concat(stats)