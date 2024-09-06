from os import name
from dash_extensions.enrich import html
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac
import feffery_utils_components as fuc
from dash import dcc, no_update, ctx, set_props, MATCH, ALL, Patch
from dash_extensions.enrich import clientside_callback, Output, Input, ClientsideFunction, callback
from dash.exceptions import PreventUpdate

import plotly.express as px
import pandas as pd
import importlib.resources as res
import anndata

import uuid
from flask_login import current_user

from stintev.utils._plot import *
from stintev.utils._methods import cellchatdb_query, bivariate_spatial_association
from stintev.components import SelectWithColor

class TabLigandReceptor():
    
    # const
    _width_sider = '15vw'
    _width_sider_collapsed = 60
    _width_drawer = 1000
    _rowHeight_plot_panel = 350
    _height_plot_panel_item: int = 300
    
    
    @staticmethod
    def error_handler(e):
        set_props(
            'notifications-container',
            {
                "children": dmc.Notification(
                    title = 'Input Error!',
                    action = 'show',
                    message = f'ligand or receptor is empty'
                )
            }
        )
    
    @staticmethod
    def empty_scatter():
        return px.scatter().update_layout(
            margin=dict(l=0, r=0, t=0, b=0),
            plot_bgcolor = '#ffffff', 
            uirevision='constant',
            legend_itemsizing = 'constant'
        ).update_xaxes(visible=False).update_yaxes(visible=False)
    
    def plot_panel(self, key: str):
        
        panel = fuc.FefferyGridItem(
            key = key,
            id={'type': 'FUCGRIDITEM-LR', 'index': key},
            className='fuc-GridItem',
            children=[
                html.Div(
                    className = 'fuc-div-plotPanel-gridItem',
                    children=[
                        dcc.Graph(
                            id={'type': 'FUCGRIDITEM_graph-LR', 'index': key},
                            figure=TabLigandReceptor.empty_scatter(),
                            responsive = True,
                            config = {
                                'autosizable': True,
                                'toImageButtonOptions': {'format': 'jpeg', 'scale': 2}
                            },
                            style = {'height': f'{self._height_plot_panel_item}px'}
                        )   
                    ]
                )
            ]
        )
        
        return panel

    def __init__(self) -> None:
        
        self.control_sample_selection = dmc.AccordionItem(
            children=[
                dmc.AccordionControl(
                    dmc.Text('Sample selection', className='dmc-Text-accordionControl'),
                ),
                dmc.AccordionPanel(
                    children = dmc.Stack(
                        gap = 0,
                        children = [
                            dmc.Text('Sample'),
                            fac.Select(
                                locale='en-us',
                                allowClear=False,
                                id = 'SELECT_sample-LR',
                                style = {'width': '100%'},
                                options=[] # to be filled by callback
                            ),
                            dmc.Text('Embedding'),
                            fac.Select(
                                locale='en-us',
                                allowClear=False,
                                id = 'SELECT_embedding-LR',
                                style = {'width': '100%'},
                                options=[] # to be filled by callback
                            )
                        ]
                    )
                )
            ],
            value='control_sample_selection'
        )
        
        self.control_pair_selection = dmc.AccordionItem(
            children=[
                dmc.AccordionControl(
                    dmc.Text('Pair selection', className='dmc-Text-accordionControl'),
                ),
                dmc.AccordionPanel(
                    children=dmc.Stack(
                        gap=5,
                        children=[
                            dmc.Text('Ligand'),
                            SelectWithColor(index='SELECTWITHCOLOR_ligand-LR').select,
                            dmc.Text('Receptor'),
                            SelectWithColor(index='SELECTWITHCOLOR_receptor-LR').select,
                            dmc.Space(h=3),
                            dmc.Button(
                                fullWidth=True,
                                children='Calculate correlation',
                                leftSection=DashIconify(icon='mdi:chart-scatter-plot', width=24),
                                id='BUTTON_calculate_correlation-LR'
                                
                            ),
                            dmc.Space(h=3),
                            dmc.Button(
                                fullWidth=True,
                                children='CellChatDB',
                                leftSection=DashIconify(icon='fluent:search-24-regular', width=24),
                                id='BUTTON_cellchatdb-LR'
                            )
                        ],
                    )
                )
            ],
            value='control_pair_selection'
        )
        
        self.drawer_cellchatdb = fac.Drawer(
            width=self._width_drawer,
            title = 'Search ligand-receptor pairs',
            id = 'DRAWER_cellchatdb-LR',
            children=dmc.Stack(gap=0, children=[
                dmc.Group(
                    grow=True,
                    children=[
                        dmc.Stack(gap=0, children=[
                            dmc.Text('Species:'),
                            fac.Select(
                                id = 'SELECT_cellchatdb_species-LR',
                                options = [
                                    {'label': i, 'value': i}
                                    for i in ['mouse', 'human', 'zebrafish']
                                ],
                                allowClear=False,
                            ),
                        ]),
                        dmc.Stack(gap=0, children=[
                            dmc.Text('Pathway:'),
                            fac.Select(
                                id = 'SELECT_cellchatdb_pathway-LR',
                                options=[], # to be filled by callback
                                allowClear=True,
                            ),
                        ]),
                    ]
                ),
                dmc.Space(h=3),
                fac.Table(
                    id='TABLE_cellchatdb-LR',
                    rowSelectionType='radio',
                    mode='server-side',
                    maxWidth=self._width_drawer-20,
                    locale='en-us',
                )
            ]),
        )
        
        self.sider = fac.Sider(
            className='fac-Sider',
            width=self._width_sider,
            theme='light',
            children=fac.Affix(
                html.Div([
                    dmc.Badge(
                        'Plot options', color='blue', variant='light', 
                        radius='xs', size='xl', fullWidth=True,
                        leftSection=DashIconify(icon='mdi:paint-outline', width=20)
                    ),
                    dmc.Accordion(
                        multiple=True,
                        chevronPosition='right',
                        variant='seperated',
                        children=[
                            self.control_sample_selection,
                            self.control_pair_selection
                        ],
                        value=['control_sample_selection', 'control_pair_selection']
                    ),
                ]),
            )
        )
        
        self.content = fac.Content(
            className='fac-Content',
            children=fuc.FefferyDiv(
                className='fuc-Grid-container-LR',
                id='DIV_grid_container-LR',
                children = fuc.FefferyGrid(
                    id='FUCGRID_content-LR',
                    className='fuc-Grid',
                    isDraggable=False,
                    children=[
                        self.plot_panel(key)
                        for key in ['ligand', 'receptor', 'correlation']
                    ],
                    layouts=[
                        dict(i=key, x=(n%3)*4, y=n//3, w=4, h=1)
                        for n,key in enumerate(['ligand', 'receptor', 'correlation'])
                    ],
                    rowHeight=self._rowHeight_plot_panel,
                    margin=[5,5],
                    containerPadding=[0,0],
                    autoSize=True,
                    isBounded=True,
                )
            )
        )
         
        self.tab = dbc.Tab(
            label = 'Ligand-receptor',
            tab_id = 'TAB-ligand_receptor',
            children=[
                fac.Layout(
                    [
                        self.sider,
                        self.content,
                    ],
                    className='fac-Layout'
                ),
                self.drawer_cellchatdb
            ]
        )
    pass

# open cellchatdb drawer
@callback(
    Output('DRAWER_cellchatdb-LR','visible'),
    Input('BUTTON_cellchatdb-LR','n_clicks'),
    prevent_initial_call=True
)
def open_drawer_cellchatdb(n_clicks):
    return True

# grid resize when window resize
clientside_callback(
    ClientsideFunction(
        namespace='LR',
        function_name='grid_container_resize',
    ),
    Output('DIV_grid_container-LR', '_height'),
    Input('DIV_grid_container-LR', '_width'),
    Input('DIV_grid_container-LR', '_height')
)

# update pathway options in drawer
@callback(
    Output('SELECT_cellchatdb_pathway-LR', 'options'),
    Input('SELECT_cellchatdb_species-LR', 'value'),
    prevent_initial_call=True
)
def update_pathway_options(species):

    path =  res.files('stintev').joinpath(
                f'resources/CellChatDB/interaction_{species}_CellChatDB.csv'
    ) # PosixPath
    df = pd.read_csv(str(path), index_col=0)
    pathways = df['pathway_name'].unique().tolist()
    
    return pathways

# update LR query result in drawer
@callback(
    Output('TABLE_cellchatdb-LR', 'columns'),
    Output('TABLE_cellchatdb-LR', 'data'),
    
    Input('SELECT_cellchatdb_species-LR', 'value'),
    Input('SELECT_cellchatdb_pathway-LR', 'value'),
    prevent_initial_call=True
)
def update_cellchatdb_table(species, pathway):

    result = cellchatdb_query(species, pathway)

    return result['columns'], result['data']

# update embedding options, LR options
@callback(
    Output('SELECT_embedding-LR', 'options'),
    Output('SELECT_embedding-LR', 'value'),
    Output({'type': 'SelectWithColor_select', 'index': 'SELECTWITHCOLOR_ligand-LR'}, 'options'),
    Output({'type': 'SelectWithColor_select', 'index': 'SELECTWITHCOLOR_ligand-LR'}, 'value'),
    Output({'type': 'SelectWithColor_select', 'index': 'SELECTWITHCOLOR_receptor-LR'}, 'options'),
    Output({'type': 'SelectWithColor_select', 'index': 'SELECTWITHCOLOR_receptor-LR'}, 'value'),
    
    Input('SELECT_sample-LR', 'value'),
)
def update_embedding_lr_options(sample):
    
    '''
        sample: str, path of selected adata
            e.g. '/data1/share/omics-viewer/stintev/datasets/public/spatial/E7.5.h5ad'
    '''
    if sample is None:
        return [], [], []
    
    adata = anndata.read_h5ad(
        sample, backed='r'
    )
    
    return (
        [ {'label': i, 'value': i} for i in adata.obsm.keys() ], None,
        [ {'label': i, 'value': i} for i in adata.var_names.to_list() ], None, 
        [ {'label': i, 'value': i} for i in adata.var_names.to_list() ], None,
    )

# update plots
@callback(
    Output({'type': 'FUCGRIDITEM_graph-LR', 'index': ALL}, 'figure'),
    
    Input('SELECT_sample-LR', 'value'),
    Input('SELECT_embedding-LR', 'value'),
    Input({'type': 'SelectWithColor_select', 'index': 'SELECTWITHCOLOR_ligand-LR'}, 'value'),
    Input({'type': 'SelectWithColor_select', 'index': 'SELECTWITHCOLOR_receptor-LR'}, 'value'),
    Input('BUTTON_calculate_correlation-LR', 'n_clicks'),
    running = [
        ( Output('BUTTON_calculate_correlation-LR', 'loading'), True, False),
        ( Output('SELECT_sample-LR', 'disabled'), True, False ),
        ( Output('SELECT_embedding-LR', 'disabled'), True, False ),
        ( Output({'type': 'SelectWithColor_select', 'index': 'SELECTWITHCOLOR_ligand-LR'}, 'disabled'), True, False ),
        ( Output({'type': 'SelectWithColor_select', 'index': 'SELECTWITHCOLOR_receptor-LR'}, 'disabled'), True, False, ) 
    ],
    prevent_initial_call=True,
    on_error = TabLigandReceptor.error_handler,
)
def update_corr_plot(sample, embedding, ligand, receptor, click_cal):

    '''
        sample: str, path of selected adata
            e.g. '/data1/share/omics-viewer/stintev/datasets/public/spatial/E7.5.h5ad'
        embedding: str, key of selected embedding
            e.g. 'spatial'
        ligand: str, selected ligand
            e.g. 'Ccl19'
        ligand: str, selected ligand
            e.g. 'Ccl19'
    '''

    if sample is None or embedding is None:
        raise PreventUpdate

    adata = anndata.read_h5ad(
        sample
    )

    tid = ctx.triggered_id

    plot_ligand = no_update
    plot_receptor = no_update
    plot_correlation = no_update

    if 'index' in tid and tid.index =='SELECTWITHCOLOR_ligand-LR':
        if ligand is not None:
            plot_ligand = plot_feature_embedding(
                adata=adata,
                preserved_cells=adata.obs_names.to_list(),
                feature = ligand,
                embedding=embedding,
            )
        else:
            plot_ligand = TabLigandReceptor.empty_scatter()

    if 'index' in tid and tid.index == 'SELECTWITHCOLOR_receptor-LR':
        if receptor is not None:
            plot_receptor = plot_feature_embedding(
                adata=adata,
                preserved_cells=adata.obs_names.to_list(),
                feature = receptor,
                embedding=embedding,
            )
        else:
            plot_receptor = TabLigandReceptor.empty_scatter()

    if tid == 'BUTTON_calculate_correlation-LR':

        if ligand is not None and receptor is not None:
            corr = bivariate_spatial_association(
                adata=adata,
                ligand=ligand,
                receptor=receptor,
                single_cell=False,
                embedding=embedding
            )
            plot_correlation = plot_feature_embedding(
                adata = adata,
                preserved_cells=adata.obs_names.to_list(),
                feature = corr,
                embedding=embedding,
                legend_title=f'{ligand}-{receptor}',
                cmap = [
                    (0.00, "#3D739C"),
                    (0.5, "#F4F4F4"),
                    (1.00, "#c85863"),
                ],
            )
        else:
            raise Exception('ligand or receptor empty')

    return [plot_ligand, plot_receptor, plot_correlation]

# sync camera between plots
@callback(
    output = {
        'graphs': Output({'type': 'FUCGRIDITEM_graph-LR', 'index': ALL}, 'figure'),
    },
    inputs = {
        'relayoutData': Input({'type': 'FUCGRIDITEM_graph-LR', 'index': ALL}, 'relayoutData'),
    },
    prevent_initial_call = True
)
def sync_camera_between_plots_LR(relayoutData):
    
    tid = ctx.triggered_id
    graph_order = ['ligand', 'receptor', 'correlation'].index(tid['index'])
    view_to_set = relayoutData[graph_order]
    patch = Patch()
    if 'scene.camera' in view_to_set:
        patch['layout']['scene']['camera'] = view_to_set['scene.camera']
    if 'scene.aspectratio' in view_to_set:
        patch['layout']['scene']['aspectmode'] = 'manual'
        patch['layout']['scene']['aspectratio'] = view_to_set['scene.aspectratio']
    
    return_graphs = [no_update]*3
    for i in range(3):
        if i != graph_order:
            return_graphs[i] = patch

    return {'graphs': return_graphs}