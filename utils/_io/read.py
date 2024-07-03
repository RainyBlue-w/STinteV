from typing import Optional, Literal
import anndata
import os

def read_dataset(
    path_dataset: str, 
    backed:  bool | Literal['r', 'r+'] | None = None,
):
    dataset = {
        file : anndata.read_h5ad(
            filename=os.path.join(path_dataset, file), 
            backed=backed
        ) 
        for file in os.listdir(path_dataset) if file.endswith('.h5ad')
    }
    return dataset

def read_txt(
    path_txt: str
):
    with open(path_txt) as f:
        content = f.read()
    return content