from dash_extensions.enrich import callback, Output, Input, State, Serverside, html, no_update, Trigger
from dash import clientside_callback, ClientsideFunction, ALL, MATCH, Patch, ctx
from dash.exceptions import PreventUpdate
import dash_mantine_components as dmc
from dash_iconify import DashIconify

import anndata
import os
import numpy as np
from typing import List, Literal
from flask_login import current_user

from stintev.config import PathConfig
from stintev.utils._plot import plot_feature_embedding, plot_metadata_embedding
from stintev.components import PlotPanel, DatasetList, PanelLinkages, DataFilter
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

    n_outputs = len(ctx.outputs_grouping['styles'])
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
def add_delete_plot_panel(add, delete, uuid_list, choosen_dataset, path_server_folder):
    
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
            DatasetList.alert_guest()
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
clientside_callback( # update PlotPanel display_idx
    ClientsideFunction(
        namespace='overview',
        function_name='update_PlotPanel_display_idx'
    ),
    Output({'type': 'PlotPanel_badge_panel_idx', 'index': ALL}, 'children'),
    Input('STORE_plotPanelsCurUUID-overview', 'data')
)

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
                                    for adata in sorted( os.listdir(dataset_dir) )  if adata.endswith('.h5ad')
                                    # if os.path.exists( os.path.join(dataset_dir, adata) )
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
    Output({'type': 'DataFilter_select_column', 'index': MATCH}, 'options'), # filter column options
    Output({'type': 'PlotPanel_item_select_embedding', 'index': MATCH}, 'options'),
    Input({'type': 'PlotPanel_item_select_info', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
    prevent_initial_call=True
)
def update_PlotPanel_column_options(info, path_sample):

    adata = anndata.read_h5ad(
        path_sample, backed='r'
    )
    
    tid = ctx.triggered_id
    
    obs_columns = adata.obs.columns.to_list()
    
    if len(tid)==1 and tid['type']=='PlotPanel_item_select_info':
        if info == 'feature':
            options_column = [{'label': column, 'value': column} for column in adata.var_names]
            return options_column, no_update, no_update
        elif info == 'metadata':
            options_column = [{'label': column, 'value': column} for column in adata.obs.columns]
            return options_column, no_update, no_update
        
    elif (len(tid)==1 and tid['type']=='PlotPanel_item_select_sample') or (len(tid)==2):
        options_embedding = [{'label': i, 'value': i} for i in list(adata.obsm.keys())]
        if info and (info=='feature'):
            options_column = [{'label': column, 'value': column} for column in adata.var_names]
            return options_column, obs_columns, options_embedding
        elif info and (info=='metadata'):
            options_column = [{'label': column, 'value': column} for column in adata.obs.columns]
            return options_column, obs_columns, options_embedding
        else:
            return no_update, no_update, options_embedding
    else:
        raise PreventUpdate

#clientside
@dashapp.callback( # update PlotPanel figure
    Output({'type': 'PlotPanel_item_graph', 'index': MATCH}, 'figure'),
    Output({'type': 'PlotPanel_store_traceNumber', 'index': MATCH}, 'data'),
    
    Input({'type': 'PlotPanel_item_select_column', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_embedding', 'index': MATCH}, 'value'),
    
    Input({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_info', 'index': MATCH}, 'value'),
    Input({'type': 'DataFilter_store_preserved_cells', 'index': MATCH}, 'data'),
    prevent_initial_call=True
)
def update_PlotPanel_figure(column, embedding, path_sample, info, preserved_cells):
    
    adata = anndata.read_h5ad(
        path_sample, backed='r'
    )
    
    if preserved_cells is None or len(preserved_cells) < 1:
        preserved_cells = adata.obs_names
    
    if info == 'feature' and column and embedding and path_sample and info:
        figure = plot_feature_embedding(adata, preserved_cells, column, embedding)
        traceNumber = 1
    elif info == 'metadata' and column and embedding and path_sample and info:
        figure = plot_metadata_embedding(adata, preserved_cells, column, embedding)
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

#region linkage

#clientside
@dashapp.callback( # update options and value for linkage panels options
    Output({'type': 'PanelLinkages_select_linkage', 'index': ALL}, 'data'),
    Output({'type': 'PanelLinkages_select_linkage', 'index': ALL}, 'value'),
    
    Input('STORE_plotPanelsCurUUID-overview', 'data'), # update when add/delete PlotPanel
    State({'type': 'PanelLinkages_select_linkage', 'index': ALL}, 'value'),

    prevent_initial_call=False
)
def update_linkage_panels_select_opitons_and_value(list_uuid, list_state_value):

    options = [{'label': f'Panel {i+1}', 'value': uid} for i,uid in enumerate(list_uuid)]
    
    values = []
    for state_value in list_state_value:
        if state_value is not None:
            values.append([ i for i in state_value if i in list_uuid])
        else:
            values.append(no_update)

    return [options]*len(ctx.outputs_list[0]), values

@dashapp.callback( # update linkage marks
    Output({'type': 'PlotPanel_linkage_marks', 'index': ALL}, 'children'),
    inputs={
        'all_inputs': {
            'list_types': Input({'type': 'PanelLinkages_select_type', 'index': ALL}, 'value'),
            'list_panels': Input({'type': 'PanelLinkages_select_linkage', 'index': ALL}, 'value'),
            'list_apply': Input({'type': 'PanelLinkages_switch_apply', 'index': ALL}, 'checked'),
        },
        'all_states': {
            'list_plotPanels_curUUID': State('STORE_plotPanelsCurUUID-overview', 'data')
        }
    },
)
def apply_linkage(all_inputs, all_states):
    
    '''ctx.args_grouping.all_inputs: 
    {
        'list_types': [
            {
                'id': {'index': '4d2cf75a435011ef959c3cecef387085', 'type': 'PanelLinkages_select_type'}, 
                'property': 'value', 'value': ['sample', 'embedding'], 
                'str_id': '{"index":"4d2cf75a435011ef959c3cecef387085","type":"PanelLinkages_select_type"}', 
                'triggered': True
            }, 
            {
                'id': {'index': '547bc360435011efab683cecef387085', 'type': 'PanelLinkages_select_type'}, 
                'property': 'value', 'value': ['column', 'view'], 
                'str_id': '{"index":"547bc360435011efab683cecef387085","type":"PanelLinkages_select_type"}', 
                'triggered': False
            }
        ], 
        'list_panels': [
            {
                'id': {'index': '4d2cf75a435011ef959c3cecef387085', 'type': 'PanelLinkages_select_linkage'}, 
                'property': 'value', 
                'value': ['4d28e2e6435011ef959c3cecef387085', '4cee8880435011ef959c3cecef387085'], 
                'str_id': '{"index":"4d2cf75a435011ef959c3cecef387085","type":"PanelLinkages_select_linkage"}', 
                'triggered': False
            }, 
            {
                'id': {'index': '547bc360435011efab683cecef387085', 'type': 'PanelLinkages_select_linkage'}, 
                'property': 'value', 
                'value': ['4d28e2e6435011ef959c3cecef387085', '4cee8880435011ef959c3cecef387085'], 
                'str_id': '{"index":"547bc360435011efab683cecef387085","type":"PanelLinkages_select_linkage"}', 
                'triggered': False
            }
        ], 
        'list_apply': [
            {
                'id': {'index': '4d2cf75a435011ef959c3cecef387085', 'type': 'PanelLinkages_switch_apply'}, 
                'property': 'checked', 'value': False, 
                'str_id': '{"index":"4d2cf75a435011ef959c3cecef387085","type":"PanelLinkages_switch_apply"}', 
                'triggered': False
            }, 
            {
                'id': {'index': '547bc360435011efab683cecef387085', 'type': 'PanelLinkages_switch_apply'}, 
                'property': 'checked', 'value': False, 
                'str_id': '{"index":"547bc360435011efab683cecef387085","type":"PanelLinkages_switch_apply"}', 
                'triggered': False
            }
        ]
    }
    '''

    color_sequence = ['red', 'blue', 'green', 'orange', 'purple', 'cyan']
    i = -1
    list_linkage_type = [item['value'] for item in ctx.args_grouping.all_inputs.list_types]
    list_selected_panels = [item['value'] for item in ctx.args_grouping.all_inputs.list_panels]
    list_linkage_apply = [item['value'] for item in ctx.args_grouping.all_inputs.list_apply]
    panels_curUUID = all_states['list_plotPanels_curUUID']
    marks = [[]] * len(panels_curUUID)
    for type, panels, apply in zip(list_linkage_type, list_selected_panels, list_linkage_apply): # 对每条linkage
        # 颜色迭代器
        i += 1
        i = i % 3
        if apply and type: # 如果应用linkage
            for panel in panels: # 对于选中的每个panel（uid）
                index = panels_curUUID.index(panel) # 找到uid对应的现有panel的位置
                print(index)
                marks[index] = marks[index] + [
                    PanelLinkages.linkage_mark(color=color_sequence[i])
                ] # 不能改成marks[index].append(...)或者marks[index]+=...，否则从marks为[[],[], ...]时每次都会会给所有panels添加mark
        else: # 不应用的话则跳过
            continue
    return marks

#clientside
@dashapp.callback( # update linked panels (sample)
    Output({'type': 'PlotPanel_item_select_sample', 'index': ALL}, 'value'),
    
    Input({'type': 'PlotPanel_item_select_sample', 'index': ALL}, 'value'),
    State('PanelLinkage_select_type', 'value'),
    State('PanelLinkage_select_linkage', 'value'),
    State('STORE_plotPanelsCurUUID-overview', 'data'),
)
def update_linked_panels_sample(
    list_selected_sample: List, 
    linkage_type: List[Literal[None, 'column','view', 'sample', 'embedding']], 
    linkage_panels: List, 
    list_cur_uuid: List
):
    
    if linkage_type is [None] or linkage_type is None:
        raise PreventUpdate

    if 'sample' in linkage_type and len(linkage_panels) >= 2:
        return_list = [no_update]*len(list_cur_uuid)
        tid = ctx.triggered_id
        sample_to_set = list_selected_sample[list_cur_uuid.index(tid['index'])]
        for i,uid in enumerate(list_cur_uuid):
            if uid in linkage_panels and uid != ctx.triggered_id['index']:
                return_list[i] = sample_to_set
        return return_list

    raise PreventUpdate

@dashapp.callback( # update linked panels (embeddings)
    Output({'type': 'PlotPanel_item_select_embedding', 'index': ALL}, 'value'),

    Input({'type': 'PlotPanel_item_select_embedding', 'index': ALL}, 'value'),
    State('PanelLinkage_select_type', 'value'),
    State('PanelLinkage_select_linkage', 'value'),
    State('STORE_plotPanelsCurUUID-overview', 'data'),
)
def update_linked_panels_embedding(
    list_selected_embedding: List, 
    linkage_type: List[Literal[None, 'column','view', 'sample', 'embedding']], 
    linkage_panels: List, 
    list_cur_uuid: List
):
    
    if linkage_type is [None] or linkage_type is None:
        raise PreventUpdate

    if 'embedding' in linkage_type and len(linkage_panels) >= 2:
        return_list = [no_update]*len(list_cur_uuid)
        tid = ctx.triggered_id
        embedding_to_set = list_selected_embedding[list_cur_uuid.index(tid['index'])]
        for i,uid in enumerate(list_cur_uuid):
            if uid in linkage_panels and uid != ctx.triggered_id['index']:
                return_list[i] = embedding_to_set
        return return_list

    raise PreventUpdate

@dashapp.callback( # update linked panels (column)
    Output({'type': 'PlotPanel_item_select_column', 'index': ALL}, 'value'),
    Output({'type': 'PlotPanel_item_select_info', 'index': ALL}, 'value'),
    
    Input({'type': 'PlotPanel_item_select_column', 'index': ALL}, 'value'),
    Input({'type': 'PlotPanel_item_select_info', 'index': ALL}, 'value'),
    State('PanelLinkage_select_type', 'value'),
    State('PanelLinkage_select_linkage', 'value'),
    State('STORE_plotPanelsCurUUID-overview', 'data'),
    prevent_initial_call=True
)
def update_linked_panels_column(list_selected_column: List, list_selected_info: List, 
                                linkage_type: List[Literal[None, 'column','view', 'sample', 'embedding']], linkage_panels: List, list_cur_uuid: List):
    
    '''
    ctx.inputs_list: [
        [
            {'id': {'index': 'c639caf83dfe11ef997e3cecef387085', 'type': 'PlotPanel_item_select_column'}, 'property': 'value', 'value': 'leiden'}, 
            {'id': {'index': 'c65229b83dfe11ef997e3cecef387085', 'type': 'PlotPanel_item_select_column'}, 'property': 'value', 'value': 'leiden'}
        ]
    ]
    list_selected_column: ['leiden', 'leiden']
    linkage_panels: ['c639caf83dfe11ef997e3cecef387085', 'c65229b83dfe11ef997e3cecef387085']
    '''
    
    if linkage_type is [None] or linkage_type is None:
        raise PreventUpdate

    if 'column' in linkage_type and len(linkage_panels) >= 2:
        return_list_column = [no_update]*len(list_cur_uuid)
        return_list_info = [no_update]*len(list_cur_uuid)
        tid = ctx.triggered_id
        column_to_set = list_selected_column[list_cur_uuid.index(tid['index'])]
        info_to_set = list_selected_info[list_cur_uuid.index(tid['index'])]
        for i,uid in enumerate(list_cur_uuid):
            if uid in linkage_panels and uid != ctx.triggered_id['index']:
                return_list_column[i] = column_to_set
                return_list_info[i] = info_to_set
        return return_list_column, return_list_info
    
    raise PreventUpdate

#clientside(set_props)
@dashapp.callback( # update linked panels (view)
    Output({'type': 'PlotPanel_item_graph', 'index': ALL}, 'figure'),
    
    Input({'type': 'PlotPanel_item_graph', 'index': ALL}, 'relayoutData'),
    State('PanelLinkage_select_type', 'value'),
    State('PanelLinkage_select_linkage', 'value'),
    State('STORE_plotPanelsCurUUID-overview', 'data'),
    prevent_initial_call=True
)
def update_linked_panels_view(
    list_relayoutData, 
    linkage_type: List[Literal[None, 'column','view', 'sample', 'embedding']], 
    linkage_panels: List, list_cur_uuid: List
):

    if linkage_type is [None] or linkage_type is None:
        raise PreventUpdate

    if 'view' in linkage_type and len(linkage_panels) >= 2:
        fig_updates = [no_update]*len(list_cur_uuid)
        tid = ctx.triggered_id
        view_to_set = list_relayoutData[list_cur_uuid.index(tid['index'])]
        for i,uid in enumerate(list_cur_uuid):
            if uid in linkage_panels and uid != ctx.triggered_id['index']:
                patch = Patch()
                if 'scene.camera' in view_to_set:
                    patch['layout']['scene']['camera'] = view_to_set['scene.camera']
                if 'scene.aspectratio' in view_to_set:
                    patch['layout']['scene']['aspectmode'] = 'manual'
                    patch['layout']['scene']['aspectratio'] = view_to_set['scene.aspectratio']
                fig_updates[i] = patch
                
        return fig_updates
                
    raise PreventUpdate

#endregion

#region filter

@dashapp.callback( # update filter column options
    Output({'type': 'DataFilter_select_column', 'index': MATCH}, 'options'),

    Input({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
)
def update_filter_column_options(path_sample):

    adata = anndata.read_h5ad(path_sample, backed='r')
    
    return adata.obs.columns.to_list()

@dashapp.callback( # update item type (numeric/categorical)
    Output({'type': 'DataFilter_select_type', 'index': MATCH}, 'value'),
    Output({'type': 'DataFilter_switch_apply', 'index': MATCH}, 'checked'),

    Input({'type': 'DataFilter_select_column', 'index': MATCH}, 'value'),
    State({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
)
def update_selectData_FilterType_3D(selected_column, path_sample):

    adata = anndata.read_h5ad(path_sample, backed='r')

    if selected_column is None:
        return no_update, False

    dtype = adata.obs[selected_column].dtype
    columnType = 'categorical' if ( dtype in [np.dtype('O'), 'category']) else 'numeric'

    return columnType, False

@dashapp.callback( # generate filter body
    Output({'type': 'DataFilter_filter_body', 'index': MATCH}, 'children'),
    
    Input({'type': 'DataFilter_select_type', 'index': MATCH}, 'value'),
    State({'type': 'DataFilter_select_column', 'index': MATCH}, 'value'),
    
    State({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
)
def dataFilter_generate_filter_body(
    type: Literal[None, 'numeric','categorical'], 
    selected_column: str,
    path_sample: str,
):
    
    adata = anndata.read_h5ad(path_sample, backed='r')
    
    if type == 'numeric':
        return DataFilter.numeric_filter(
            index=ctx.triggered_id['index'], 
            column=selected_column, 
            min = adata.obs[selected_column].min(), 
            max = adata.obs[selected_column].max()
        )
    
    if type == 'categorical':
        return DataFilter.categorical_filter(
            index=ctx.triggered_id['index'], 
            options=sorted(adata.obs[selected_column].unique())
        )
    
    raise PreventUpdate

#clientside
@dashapp.callback( # update STORE(serverside) for preserved_cells
    Output({'type': 'DataFilter_store_preserved_cells', 'index': MATCH}, 'data'),
    Output({'type': 'DataFilter_text_number', 'index': MATCH}, 'children'),
    
    Input({'type': 'DataFilter_switch_apply', 'index': MATCH}, 'checked'),
    State({'type': 'DataFilter_select_column', 'index': MATCH}, 'value'),
    State({'type': 'DataFilter_select_type', 'index': MATCH}, 'value'),
    Input({'type': 'DataFilter_numberInput_left', 'index': MATCH}, 'value'),
    Input({'type': 'DataFilter_numberInput_right', 'index': MATCH}, 'value'),
    Input({'type': 'DataFilter_transfer', 'index': MATCH}, 'targetKeys'),
    Input({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
)
def dataFilter_apply_filter_store_preserved_cells(
    checked: bool,
    selected_column: str,
    type: Literal[None, 'numeric','categorical'],
    min: float,
    max: float,
    transfer: str,
    path_sample: str,
):
    
    adata = anndata.read_h5ad(path_sample, backed='r')
    
    if checked:
        if type == 'numeric':
            preserved_cells = adata.obs.index[
                (adata.obs[selected_column] >= min) & (adata.obs[selected_column] <= max)
            ].tolist()
            return Serverside(preserved_cells), f'Selected cells: {len(preserved_cells)}'
        
        if type == 'categorical':
            preserved_cells = adata.obs.index[
                adata.obs[selected_column].isin(transfer)
            ].tolist()
            return Serverside(preserved_cells), f'Selected cells: {len(preserved_cells)}'
    else:
        return Serverside(adata.obs_names.to_list()), 'Selected cells: 0'

# endregion