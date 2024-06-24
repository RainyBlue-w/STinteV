
import anndata
import pandas as pd

import plotly.express as px
import plotly.graph_objects as go

from typing import List, Dict
from pydantic import BaseModel

# class for the verification of parameter passing
class ParamsPlotFeatureEmbedding(BaseModel):
    feature: str
    embedding: str
    sort: bool = True
    ascending: bool = True
    cmap: List = [
        (0.00, "#F4F4F4"),
        (0.05, "#F4F4F4"),
        (1.00, "#225EA8")
    ]

class ParamsPlotMetadataEmbedding(BaseModel):
    column: str
    embedding: str
    cmap: Dict | None = None

# plot functions

def plot_feature_embedding(
    adata: anndata.AnnData,
    feature: str,
    embedding: str,
    sort: bool = True,
    ascending: bool = True,
    cmap: List = [
        (0.00, "#F4F4F4"),
        (0.05, "#F4F4F4"),
        (1.00, "#225EA8")
    ],
    marker_size: float = 2.0,
    **kws
) -> go.Figure :
    
    if adata.obsm[embedding].shape[1] == 2:
        plot = _plot_feature_embedding_2d(adata, feature, embedding, sort, ascending, cmap, marker_size, **kws)
    elif adata.obsm[embedding].shape[1] == 3:
        plot = _plot_feature_embedding_3d(adata, feature, embedding, sort, ascending, cmap, marker_size, **kws)
    else:
        raise ValueError(f"The embedding '{embedding}' seems to be neither 2D nor 3D")
    
    return plot

def plot_metadata_embedding(
    adata: anndata.AnnData,
    column: str,
    embedding: str,
    cmap: Dict = None,
    marker_size: float = 2.0,
    color_discrete_sequence: List = px.colors.qualitative.Alphabet,
    **kws
) -> go.Figure:
    
    if adata.obsm[embedding].shape[1] == 2:
        plot = _plot_metadata_embedding_2d(adata, column, embedding, cmap, marker_size, color_discrete_sequence=color_discrete_sequence, **kws)
    elif adata.obsm[embedding].shape[1] == 3:
        plot = _plot_metadata_embedding_3d(adata, column, embedding, cmap, marker_size, color_discrete_sequence=color_discrete_sequence, **kws)
    else:
        raise ValueError(f"The embedding '{embedding}' seems to be neither 2D nor 3D")
    
    return plot

def _plot_feature_embedding_2d(
    adata: anndata.AnnData,
    feature: str,
    embedding: str,
    sort: bool = True,
    ascending: bool = True,
    cmap: List = [
        (0.00, "#F4F4F4"),
        (0.05, "#F4F4F4"),
        (1.00, "#225EA8")
    ],
    marker_size: float = 2.0,
    **kws
) -> go.Figure :
    
    """
    plot feature on 2D-embedding
    """

    embedding = pd.DataFrame(adata.obsm[embedding], index=adata.obs_names, columns=['X','Y'])
    pdf = pd.concat([embedding, adata[:,feature].to_df()], axis=1)
    if sort is True:
        pdf = pdf.sort_values(by=feature, ascending=ascending)
    plot = px.scatter(
        data_frame = pdf,
        x = 'X', y = 'Y', color = feature,
        color_continuous_scale = cmap, render_mode='webgl',
        **kws
    )
    plot.update_yaxes(visible=False).update_xaxes(visible=False)
    plot.update_traces(marker_size=marker_size, marker_opacity=1)
    plot.update_layout(
        margin=dict(l=0, t=0, b=0),
        plot_bgcolor = '#ffffff', 
        uirevision='constant',
        legend_itemsizing = 'constant',
        legend=dict(
            title='',
            orientation='v',
            yanchor='middle',
            xanchor='right',  # 设置图例在右侧
            y=0.5,
            x=1,  # 调整图例在横向的位置
        ),
    )
    plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],font_size = 20))
    
    return plot

def _plot_feature_embedding_3d(
    adata: anndata.AnnData,
    feature: str,
    embedding: str,
    sort: bool = True,
    ascending: bool = True,
    cmap: List = [
        (0.00, "#F4F4F4"),
        (0.05, "#F4F4F4"),
        (1.00, "#225EA8")
    ],
    marker_size: float = 2.0,
    **kws
) -> go.Figure :

    """
    plot feature on 3D-embedding
    """
    pdf = pd.DataFrame(adata.obsm[embedding], index=adata.obs_names, columns=['X','Y','Z'])
    pdf = pd.concat([pdf, adata[:,feature].to_df()], axis=1)
    if sort is True:
      pdf = pdf.sort_values(by=feature, ascending=ascending)
    plot = px.scatter_3d(
        data_frame = pdf,
        x = 'X', y = 'Y', z = 'Z',
        color = feature,
        color_continuous_scale = cmap,
        **kws
    )
    plot.update_traces(marker_size=marker_size, marker_opacity=1)
    plot.update_layout(
      margin=dict(l=0, r=0, t=0, b=0),
      plot_bgcolor = '#ffffff', 
      uirevision='constant',
      coloraxis = {
        'colorbar' : {'tickformat': '4.2f'}
      },
      scene = dict(
          xaxis = dict(backgroundcolor='white', showbackground=True, zerolinecolor='gray', autorange=True,
                       gridcolor='gray', nticks=6),
          yaxis = dict(backgroundcolor='white', showbackground=True, zerolinecolor='gray', autorange=True,
                       gridcolor='gray', nticks=6),
          zaxis = dict(backgroundcolor='white', showbackground=True, zerolinecolor='gray', autorange=True,
                       gridcolor='gray', nticks=6),
          bgcolor = 'white',
          camera = dict(projection = dict(type='orthographic') ),
          aspectmode = 'data'
      )
    )
    # plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],
    #                                            font_size = 20)) 
    return plot

def _plot_metadata_embedding_2d(
    adata: anndata.AnnData,
    column: str,
    embedding: str,
    cmap: Dict | None = None,
    marker_size: float = 2.0,
    **kws
) -> go.Figure:
    
    pdf = pd.DataFrame(adata.obsm[embedding], index=adata.obs_names, columns=['X', 'Y'])
    pdf = pd.concat([ pdf, adata.obs[column] ], axis=1)
    pdf = pdf.sort_values(by=column)
    plot = px.scatter(
        data_frame = pdf,
        x = 'X', y = 'Y', color = column,
        color_discrete_map = cmap,
        **kws
    )
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=marker_size, marker_opacity=1)
    plot.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        plot_bgcolor = '#ffffff', 
        uirevision='constant',
        legend_itemsizing = 'constant'
    )
    plot.for_each_annotation(
        lambda a: a.update(
            text=a.text.split("=")[-1],font_size = 20
        )
    )

    return plot

def _plot_metadata_embedding_3d(
    adata: anndata.AnnData,
    column: str,
    embedding: str,
    cmap: Dict | None = None,
    marker_size: float = 2.0,
    **kws
) -> go.Figure:
    
    pdf = pd.DataFrame(adata.obsm[embedding], index=adata.obs_names, columns=['X', 'Y', 'Z'])
    pdf = pd.concat([ pdf, adata.obs[column] ], axis=1)
    pdf = pdf.sort_values(by=column)
    plot = px.scatter_3d(
        data_frame = pdf,
        x = 'X', y = 'Y', z='Z', color = column,
        color_discrete_map = cmap,
        **kws
    )
    plot.update_traces(marker_size=marker_size, marker_opacity=1)
    plot.update_layout(
      margin=dict(l=0, r=0, t=0, b=0),
      plot_bgcolor = '#ffffff', 
      uirevision='constant',
      legend_itemsizing = 'constant',
      coloraxis = {
        'colorbar' : {'tickformat': '4.2f'}
      },
      scene = dict(
          xaxis = dict(backgroundcolor='white', showbackground=True, zerolinecolor='gray', autorange=True,
                       gridcolor='gray', nticks=6),
          yaxis = dict(backgroundcolor='white', showbackground=True, zerolinecolor='gray', autorange=True,
                       gridcolor='gray', nticks=6),
          zaxis = dict(backgroundcolor='white', showbackground=True, zerolinecolor='gray', autorange=True,
                       gridcolor='gray', nticks=6),
          bgcolor = 'white',
          camera = dict(projection = dict(type='orthographic') ),
          aspectmode = 'data'
      )
    )

    return plot


if __name__ == '__main__':
    
    from sys import getsizeof
    
    adata = anndata.read_h5ad('/rad/share/omics-viewer/spatial/matrix_data/embryo_3-2-E8.0_min400_Ann_HC0.5.h5ad', backed='r')
    
    _plot_feature_embedding_2d(adata, feature='T', embedding='X_sagittal')

    