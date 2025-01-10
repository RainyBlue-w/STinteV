import spatialdm as sdm
import anndata
from typing import List, Literal, Dict
import pandas as pd
import importlib.resources as res

def bivariate_spatial_association(
    adata_path: str,
    ligand: str,
    receptor: str,
    embedding: str = 'spatial',
    l: float = 1.2,
    cutoff: float = 0.2,
    n_neighbors: int | None = None,
    n_nearest_neighbors: int = 6,
    single_cell: bool = False,
    method: Literal['z-score', 'permutation'] = 'z-score',
) -> List:
    
    adata = anndata.read_h5ad(
        adata_path, backed=False
    )
    sdm.weight_matrix(
        adata, l=l, embedding=embedding, cutoff=cutoff, single_cell=single_cell,
        n_neighbors=n_neighbors, n_nearest_neighbors=n_nearest_neighbors
    )
    adata.uns['ligand'] = pd.DataFrame(
        data = {
            'Ligand0': ligand,
            'Ligand1': None,
            'Ligand2': None,
        },
        index = [f'{ligand}_{receptor}']
    )
    adata.uns['receptor'] = pd.DataFrame(
        data = {
            'Receptor0': receptor,
            'Receptor1': None,
            'Receptor2': None,
        },
        index = [f'{ligand}_{receptor}']
    )
    adata.uns['mean'] = 'algebra'
    adata.uns['geneInter'] = pd.DataFrame(
        data = {
            'interaction_name': f'{ligand}_{receptor}',
            'annotation': 'pseudo_interaction'
        },
        index = [f'{ligand}_{receptor}']
    )
    
    sdm.spatialdm_local(
        adata, n_perm=1000, nproc=4,
        specified_ind=[f'{ligand}_{receptor}'],
        method=method
    )
    sdm.sig_spots(adata, method=method, fdr=False, threshold=0.05)
    
    if method == 'z-score':
        return ( 1 - adata.uns['local_z_p'].loc[f'{ligand}_{receptor}'] ).to_list()
    if method == 'permutation':
        return ( 1 - adata.uns['local_perm_p'].loc[f'{ligand}_{receptor}'] ).to_list()

def cellchatdb_query(
    species: Literal['human', 'mouse', 'zebrafish'],
    pathway: str | None = None,
    mode: Literal['local', 'online'] = 'local',
) -> Dict:
    
    if mode == 'local':
        geneInter =  pd.read_csv(
            str( res.files('stintev').joinpath(
                f'resources/CellChatDB/interaction_{species}_CellChatDB.csv'
            ) ), 
            index_col=0
        )
    else:
        if species == 'mouse':
            geneInter = pd.read_csv('https://figshare.com/ndownloader/files/36638919', index_col=0)
        elif species == 'human':
            geneInter = pd.read_csv('https://figshare.com/ndownloader/files/36638943', header=0, index_col=0)
        elif species == 'zebrafish':
            geneInter = pd.read_csv('https://figshare.com/ndownloader/files/38756022', header=0, index_col=0)

    if pathway is not None:
        geneInter = geneInter.query('pathway_name == @pathway')
    
    columns = [
        {'title': i, 'dataIndex': i}
        for i in geneInter.columns.to_list()
    ]
    data = [
        {
            **record,
            'key': record['interaction_name']
        }
        for i,record in enumerate(geneInter.to_dict('records'))
    ]

    return {
        'columns': columns, 
        'data': data, 
        # 'filterOptions': {
        #     column['dataIndex']: {'filterMode': 'keyword'} 
        #     for column in columns
        # }
    }