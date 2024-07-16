import subprocess
import sys
import os
import shutil
import pandas as pd
from scipy.io import mmread
import anndata
from scipy import sparse


def rds_to_h5ad(
        file_path: str
):
    output_file_path = os.path.splitext(file_path)[0] + ".h5ad"
    file_path_tmp = file_path+"_tmp/"
    reduction_path = file_path_tmp+"reduction/"
    counts_path = file_path_tmp+"counts.mtx"
    gene_path = file_path_tmp+"gene.csv"
    meta_path = file_path_tmp+"meta.csv"
    if run_r_script(file_path):
        meta = pd.read_csv(meta_path)
        gene = pd.read_csv(gene_path)
        reductions = {}
        for file in os.listdir(reduction_path):
            key = file.split(".")[0]
            reductions[key] = pd.read_csv(reduction_path+file)
            reductions[key].set_index(reductions[key].columns[0], inplace=True)
        counts = mmread(counts_path)
        meta.set_index(meta.columns[0], inplace=True)
        gene.set_index(gene.columns[0], inplace=True)
        adata = anndata.AnnData(X=counts.transpose(), obs=meta, var=gene)
        for key in reductions:
            adata.obsm[key] = reductions[key].values
        if sparse.isspmatrix_coo(adata.X):
            adata.X = adata.X.tocsr()
        if '_index' in adata.var.columns:
            adata.var.rename(columns={'_index': 'index'}, inplace=True)
        adata.write_h5ad(output_file_path, compression='gzip')
        
    # 删除原始rds以及构建h5ad产生的临时文件
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
        except OSError as e:
            print(f"Error: {e}")

    if os.path.exists(file_path_tmp):
        try:
            shutil.rmtree(file_path_tmp)
        except OSError as e:
            print(f"Error: {e}")

def run_r_script(
        file: str
):
    try:
        result = subprocess.run(
            ["Rscript", "stintev/utils/_io/parseRDS.r", file],
            check=True,
            capture_output=True,
            text=True
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"\033[91mError: {e.stderr}\033[0m", file=sys.stderr)
        return False