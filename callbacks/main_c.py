from dash_extensions.enrich import Output, Input, State, Serverside, html, no_update
from dash_extensions.enrich import callback, clientside_callback, ClientsideFunction, callback_context
from dash import ALL, MATCH, Patch, ctx
from dash.exceptions import PreventUpdate
import dash_mantine_components as dmc
from dash_iconify import DashIconify

import anndata
import os
from typing import Literal
from flask_login import current_user

from stintev.config import PathConfig
from stintev.utils._plot import plot_feature_embedding, plot_metadata_embedding
from stintev.components import PlotPanel, DatasetList
from stintev.server import dashapp

#region update overview-grid layout

@dashapp.callback( # update the rowHeight for PlotPanel-items
    output = dict(
        styles = Output({'type': 'PlotPanel_item_graph', 'index': ALL}, 'style'),
        rowHeight = Output('FUCGRID_content-overview', 'rowHeight')
    ),
    inputs = [
        Input('NUMBERINPUT_setting_panels-overview', 'value')
    ],
    prevent_initial_call = False,
    suppress_callback_exceptions=True
)
def update_height_for_plot_panel_items(height):

    n_outputs = len(callback_context.outputs_grouping['styles'])
    return {
        'styles': [ {'height': [f'{height-75}px']} ] * n_outputs,
        'rowHeight': height
    }

@dashapp.callback( # add & delete PlotPanel
    Output('FUCGRID_content-overview', 'children', allow_duplicate=True),
    Output('STORE_plotPanelsCurUUID-overview', 'data', allow_duplicate=True),
    Output('FUCGRID_content-overview', 'layouts'),
    Output('BUTTON_setting_panels_add-overview', 'disabled'),
   
    Input('BUTTON_setting_panels_add-overview', 'n_clicks'),
    Input({'type': 'PlotPanel_item_button_delete', 'index': ALL}, 'n_clicks'),
    State('STORE_plotPanelsCurUUID-overview', 'data'),
    
    State('STORE_choosen_dataset-dataset', 'data'),
    State('STORE_server_folder-dataset', 'data'),
   prevent_initial_call=True
)
def add_plot_panel(add, delete, uuid_list, choosen_dataset, path_server_folder):
    
    import uuid
    
    tid = ctx.triggered_id
    children_grid = Patch()
    
    if tid == 'BUTTON_setting_panels_add-overview':
        if add and (len(uuid_list) <= 5):
            options = []
            for i in choosen_dataset:
                if len(i['choosen']) > 0 : 
                    for dataset in i['choosen']:
                        options.append(
                            {
                                'group': f'{dataset}-{i["group"]}',
                                'options': [
                                    {
                                        'label': adata,
                                        'value': os.path.join(path_server_folder, 'datasets', i['group'], dataset, adata)
                                    }
                                    for adata in os.listdir(
                                        os.path.join(path_server_folder, 'datasets', i['group'], dataset)
                                    )
                                ]
                            }
                        )
            
            next_PlotPanel = PlotPanel(
                uuid.uuid1().hex,
                init_samples=options,
            )
            children_grid.append(
                next_PlotPanel.grid_item
            )
            uuid_list.append(
                next_PlotPanel._index
            )
            
    elif ('type' in tid) and (tid['type'] == 'PlotPanel_item_button_delete') and (len(uuid_list) > 1):
        del_index = tid['index']
        for i, index in enumerate(uuid_list):
            if index == del_index:
                del children_grid[i]
                del uuid_list[i]
    else:
        raise PreventUpdate

    layouts = []
    for i, index in enumerate(uuid_list):
        layouts.append(
            dict(
                i=index, x = (i%3)*16, y = i//3, 
                w=16, h=1, maxH=1
            ) 
        )

    if len(uuid_list) == 6:
        ban_button = True
    else:
        ban_button = False

    return  children_grid, uuid_list, layouts, ban_button

#endregion

#region generate datalist for tab-private

@dashapp.callback(
    Output('TABS_panel_private-dataset', 'children'),
    Input('url', 'pathname'),
    prevent_initial_call=False,
)
def generate_datasetlist_private(url):
    if url == '/':
        return (
            DatasetList(
                path_data_folder=os.path.join(PathConfig.DATA_PATH,'datasets','private', current_user.username),
                group='private'
            ).list
            if current_user.is_authenticated else 
            DatasetList.alert_guest()
        )
    else:
        raise PreventUpdate
    
@dashapp.callback( # dataset创建完成后刷新private datasetlist
    Output('TABS_panel_private-dataset', 'children'),
    Input('STORE_dataset_name_upload-dataset', 'data'),
    
    prevent_initial_call=True
)
def create_refresh_datasetlist_private(dataset_name):
    if dataset_name:
        return (
            DatasetList(
                path_data_folder=os.path.join(PathConfig.DATA_PATH, 'datasets','private', current_user.username),
                group='private'
            ).list
            if current_user.is_authenticated else 
            DatasetList.alert()
        )
    else:
        raise PreventUpdate
    
@dashapp.callback( # 文件上传后刷新private datasetlit
    Output('TABS_panel_private-dataset', 'children'),
    Input('UPLOAD_dataset-dataset', 'lastUploadTaskRecord'),
)
def upload_refresh_datasetlist_private(upload):
    if upload:
        return [
            html.Div(
                id = f'TabsTabDataSet-contetn-private',
                children=[
                    DatasetList(
                        path_data_folder=os.path.join(PathConfig.DATA_PATH,'datasets','private', current_user.username),
                        group='private'
                    ).list,
                ]
            )
        ]
    else:
        raise PreventUpdate

#endregion

#region update PlotPanel

@dashapp.callback( # update dataset for tab_overview
    Output('STORE_choosen_dataset-dataset', 'data'),
    Input({'type': 'DatasetList-checkboxGroup', 'index': ALL}, 'value'),    
    prevent_initial_call=False
)
def update_choosen_dataset_for_tab_overview(value):
    
    inputs_list = ctx.inputs_list
    
    '''
    inputs_list: [
        [
            {
                'id': 
                    {'index': 'public', 'type': 'DatasetList-checkboxGroup'}, 
                'property': 'value', 
                'value': ['spatial']
            }, 
            {
                'id': {'index': 'private', 'type': 'DatasetList-checkboxGroup'}, 
                'property': 'value',
                'value': ['test2']
            }
        ]
    ]
    '''
    
    store = [ 
        dict(
            group=i['id']['index'], choosen=i['value']
        ) 
        for i in inputs_list[0] if 'value' in i
    ]
    
    '''
    store: [
        {'group': 'public', 'choosen': ['spatial']},
        {'group': 'private', 'choosen': ['test2']}
    ]
    '''
    
    return store

@dashapp.callback( # load choosen datasets in tab_dataset

    Output({'type': 'PlotPanel_item_select_sample', 'index': ALL}, 'options'),
    
    Input('STORE_choosen_dataset-dataset', 'data'),
    State('STORE_server_folder-dataset', 'data'),
    prevent_initial_call=False
)
def load_choosen_datasets(choosen_dataset, path_server_folder):
    
    ''' choosen_dataset: [
        {'group': 'public', 'choosen': ['spatial']},
        {'group': 'private', 'choosen': ['test2']}
    ]
    '''
    
    ''' ctx.outputs_list: 
    [
        {
            'id': {'index': '6579dbb22fb811ef84d03cecef387085', 'type': 'PlotPanel_item_select_sample'}, 
            'property': 'data'
        }, 
        {
            'id': {'index': '663b90b82fb811ef84d03cecef387085', 'type': 'PlotPanel_item_select_sample'}, 
            'property': 'data'
        }
    ]
    '''
    
    ''' output: options = [
        {
            'group': 'spatial',
            'items': [
                {'label': 'E7.5.h5ad', 'value': '/rad/share/omics-viewer/stintev/datasets/public/spatial/E7.5.h5ad'}, 
                {'label': E7.75.h5ad', 'value': '/rad/share/omics-viewer/stintev/datasets/public/spatial/E7.75.h5ad'},
                {'label': E8.0.h5ad', 'value': '/rad/share/omics-viewer/stintev/datasets/public/spatial/E8.0.h5ad'},
            ]
        },
        {...}
    ]
    '''

    options = []
    for i in choosen_dataset:
        if len(i['choosen']) > 0: 
            for dataset in i['choosen']:
                if  i['group'] == 'public':
                    dataset_dir = os.path.join(path_server_folder, 'datasets', i['group'], dataset)
                    if os.path.exists(dataset_dir):
                        options.append(
                            {
                                'group': f'{dataset}-{i["group"]}',
                                'options': [
                                    {
                                        'label': adata,
                                        'value': os.path.join(dataset_dir, adata)
                                    }
                                    for adata in sorted( os.listdir(dataset_dir) ) 
                                    if os.path.exists( os.path.join(dataset_dir, adata) )
                                ]
                            }
                        )
                elif i['group'] == 'private':
                    dataset_dir = os.path.join(path_server_folder,'datasets',i['group'],current_user.username,dataset)
                    if os.path.exists(dataset_dir):
                        options.append(
                            {
                                'group': f'{dataset}-{i["group"]}',
                                'options': [
                                    {
                                        'label': adata,
                                        'value': os.path.join(dataset_dir, adata)
                                    }
                                    for adata in sorted( os.listdir(dataset_dir) ) 
                                    if os.path.exists( os.path.join(dataset_dir, adata) )
                                ]
                            }
                        )

    all_options = [ options for i in range(len(ctx.outputs_list))]
    
    return all_options

@dashapp.callback( # PlotPanel updates options when sample/info update
    Output({'type': 'PlotPanel_item_select_column', 'index': MATCH}, 'options'),
    Output({'type': 'PlotPanel_item_select_embedding', 'index': MATCH}, 'options'),
    Input({'type': 'PlotPanel_item_select_info', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
    prevent_initial_call=True
)
def update_PlotPael_column_options(info, path_sample):

    adata = anndata.read_h5ad(
        path_sample, backed='r'
    )
    
    tid = ctx.triggered_id
    
    
    if len(tid)==1 and tid['type']=='PlotPanel_item_select_info':
        if info == 'feature':
            options_column = [{'label': column, 'value': column} for column in adata.var_names]
            return options_column, no_update
        elif info == 'metadata':
            options_column = [{'label': column, 'value': column} for column in adata.obs.columns]
            return options_column, no_update
        
    elif (len(tid)==1 and tid['type']=='PlotPanel_item_select_sample') or (len(tid)==2):
        options_embedding = [{'label': i, 'value': i} for i in list(adata.obsm.keys())]
        if info and (info=='feature'):
            options_column = [{'label': column, 'value': column} for column in adata.var_names]
            return options_column, options_embedding
        elif info and (info=='metadata'):
            options_column = [{'label': column, 'value': column} for column in adata.obs.columns]
            return options_column, options_embedding
        else:
            return no_update, options_embedding
    else:
        raise PreventUpdate

@dashapp.callback( # update PlotPanel figure
    Output({'type': 'PlotPanel_item_graph', 'index': MATCH}, 'figure'),
    Output({'type': 'PlotPanel_store_traceNumber', 'index': MATCH}, 'data'),
    
    Input({'type': 'PlotPanel_item_select_column', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_embedding', 'index': MATCH}, 'value'),
    
    Input({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_info', 'index': MATCH}, 'value'),
    prevent_initial_call=True
)
def update_PlotPanel_figure(column, embedding, path_sample, info):
    
    adata = anndata.read_h5ad(
        path_sample, backed='r'
    )

    if info == 'feature' and column and embedding and path_sample and info:
        figure = plot_feature_embedding(adata, column, embedding)
        traceNumber = 1
    elif info == 'metadata' and column and embedding and path_sample and info:
        figure = plot_metadata_embedding(adata, column, embedding)
        traceNumber = len(adata.obs[column].unique())
    else:
        raise PreventUpdate

    return figure, traceNumber

@dashapp.callback( # update marker_size panel
    Output({'type': 'PlotPanel_item_graph', 'index': MATCH}, 'figure'),
    Input({'type': 'PlotPanel_item_pointSize', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_store_traceNumber', 'index': MATCH}, 'data')
)
def update_PlotPanel_pointSize_panel(pointSize, traceNumber):
    patch = Patch()
    for i in range(traceNumber):
        patch['data'][i]['marker']['size'] = pointSize        
    return patch

@dashapp.callback( # update marker_size global
    Output({'type': 'PlotPanel_item_pointSize', 'index': ALL}, 'value'),
    Input('NUMBERINPUT_scatter3dPointsize_3D', 'value'),
)
def update_PlotPanel_pointSize_global(pointSize):
    return [pointSize]*len(ctx.outputs_list)

#endregion

