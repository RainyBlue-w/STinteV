import anndata
import os

def read_dataset(path_dataset: str):
    dataset = [
        anndata.read_h5ad(
            filename=os.path.join(path_dataset, file), 
            backed='r'
        ) 
        for file in os.listdir(path_dataset) if file.endswith('.h5ad')
    ]
    return dataset

