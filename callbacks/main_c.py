from dash_extensions.enrich import callback, Output, Input, State, Serverside, html, no_update, Trigger
from dash import clientside_callback, ClientsideFunction, ALL, MATCH, Patch, ctx
from dash.exceptions import PreventUpdate
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash

import anndata
import os
from fastapi import background
import numpy as np
from typing import List, Literal, Dict
from flask_login import current_user

from stintev.config import PathConfig
from stintev.utils._plot import plot_feature_embedding, plot_metadata_embedding
from stintev.components import PlotPanel, DatasetList, PanelLinkages, DataFilter
from stintev.server import dashapp
from stintev.utils._io import user_rds_path, parse_rds, convert_to_h5ad, delete_relate_tmpFiles
from stintev.utils._plot import assign_colors

#region update overview-grid layout

@dashapp.callback( # update the rowHeight for PlotPanel-items
    output = dict(
        style_graph = Output({'type': 'PlotPanel_item_graph', 'index': ALL}, 'style'),
        style_legend = Output({'type': 'PlotPanel_item_categoriesLegend', 'index': ALL}, 'style'),
        rowHeight = Output('FUCGRID_content-overview', 'rowHeight')
    ),
    inputs = [
        Input('NUMBERINPUT_setting_panels-overview', 'value')
    ],
    prevent_initial_call = False,
    suppress_callback_exceptions=True
)
def update_height_for_plot_panel_items(height):

    n_outputs = len(ctx.outputs_grouping['style_graph'])
    return_list =  {
        'style_graph' : [ {'height': [f'{height-75}px']} ] * n_outputs,
        'style_legend' : [ {'height': [f'{height-75}px'], 'overflow-y': 'auto'} ] * n_outputs,
        'rowHeight' : height,
    }
    return return_list

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
                path_data_folder=os.path.join(PathConfig.DATA_PATH,'datasets','private', current_user.id),
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
                path_data_folder=os.path.join(PathConfig.DATA_PATH, 'datasets','private', current_user.id),
                group='private'
            ).list
            if current_user.is_authenticated else 
            DatasetList.alert_guest()
        )
    else:
        raise PreventUpdate    

# 关闭转换窗口后清理临时文件
@dashapp.callback(
    Input('DIALOG_convert_rds', 'opened'),
    prevent_initial_call=True
) 
def clear_tmp_files(opened):
    if not opened:
        delete_relate_tmpFiles(user_rds_path[current_user.id])
    
# 监控step3状态变化，激活时开始转换
@dashapp.callback(
    Output('TABS_panel_private-dataset', 'children'),
    Output('STEPPER_convert_rds', 'active'),
    Output('STEP_convert_rds', 'loading'),
    Input('STEP_convert_rds', 'loading'),
    State('CHECKBOX_checked_metadata', 'value'),
    prevent_initial_call=True
) 
def start_convert_rds_to_h5ad(loading, value):
    if loading:
        path_folder = os.path.join(PathConfig.DATA_PATH,'datasets','private', current_user.id)
        convert_to_h5ad(user_rds_path[current_user.id], value)
        return [
            html.Div(
                id = f'TabsTabDataSet-contetn-private',
                children=[
                    DatasetList(
                        path_data_folder=path_folder,
                        group='private'
                    ).list,
                ]
            )
        ], 3, False
    else:
        raise PreventUpdate

# 提交元数据后跳转转换状态
@dashapp.callback(
    Output('STEPPER_convert_rds', 'active'),
    Output('STEP_convert_rds', 'loading'),
    Input('BUTTON_submit_metadata', 'n_clicks'),
    prevent_initial_call=True
)
def skip_to_rds_convert(n_clicks):
    if n_clicks is None:
        raise PreventUpdate
    else:
        return 2, True

# 监控stepper变化，启动rds文件解析工作
@dashapp.callback(
    Output('STEPPER_convert_rds', 'active'),
    Output('STEP_parse_rds', 'loading'),
    Output('STEP_parse_rds', 'color'),
    Output('STEP_parse_rds_result', 'children'),
    Output('CHECKBOX_checked_metadata', 'children'),
    Input('STEP_parse_rds', 'loading'),
    prevent_initial_call=True
)
def parse_rds_file(loading):
    if loading:
        success, metadata = parse_rds(user_rds_path[current_user.id])
        if success:
            children=dmc.Group(
                [
                    dmc.Checkbox(label=meta, value=meta) for meta in metadata
                ],
                mt=10,
            ),
            return 1, False, dash.no_update, dash.no_update, children
        return 0, False, 'red', metadata, dash.no_update
    else:
        raise PreventUpdate

# 根据上传文件类型决定是否弹出文件类型转换对话框 
@dashapp.callback( # 文件上传后刷新private datasetlit
    Output('TABS_panel_private-dataset', 'children'),
    Output('DIALOG_convert_rds', 'opened'),
    Output('STEP_parse_rds', 'loading'),
    Output('STEP_parse_rds', 'color'),
    Output('STEP_parse_rds_result', 'children'),
    Output('STEPPER_convert_rds', 'activate'),
    Input('UPLOAD_dataset-dataset', 'lastUploadTaskRecord'),
)
def upload_refresh_datasetlist_private(upload):
    if upload:
        path_folder = os.path.join(PathConfig.DATA_PATH,'datasets','private', current_user.id)
        if upload['fileName'].endswith(".rds"):
            rds_file_path = os.path.join(path_folder, upload['taskId'], upload['fileName'])
            user_rds_path[current_user.id] = rds_file_path
            return dash.no_update, True, True, 'blue', '', 0
        return [
            html.Div(
                id = f'TabsTabDataSet-contetn-private',
                children=[
                    DatasetList(
                        path_data_folder=path_folder,
                        group='private'
                    ).list,
                ]
            )
        ], False, False, dash.no_update, dash.no_update, 0
    else:
        raise PreventUpdate

#endregion

#region choose dataset

@dashapp.callback( # load choosen datasets
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

@dashapp.callback( # updates sample select options
    output = {
        'overview_panels': Output({'type': 'PlotPanel_item_select_sample', 'index': ALL}, 'options'),
        'LR_panels': Output('SELECT_sample-LR', 'options'),
    },
    inputs = {
        'choosen_dataset': Input('STORE_choosen_dataset-dataset', 'data'),
        'path_server_folder': State('STORE_server_folder-dataset', 'data'),
    },
    prevent_initial_call = True
)
def load_choosen_datasets(
    choosen_dataset: List[dict], 
    path_server_folder: str
):
    
    ''' choosen_dataset: [
        {'group': 'public', 'choosen': ['spatial']},
        {'group': 'private', 'choosen': ['test2']}
    ]
    '''
    
    ''' ctx.outputs_grouping: {
        'overview_panels': [
            {'id': {'index': '639aeabcb87e11efa5a13cecef387084', 'type': 'PlotPanel_item_select_sample'}, 'property': 'options'}, 
            {'id': {'index': '63b1b97cb87e11efa5a13cecef387084', 'type': 'PlotPanel_item_select_sample'}, 'property': 'options'}
        ], 
        'LR_panels': {'id': 'SELECT_sample-LR', 'property': 'options'}}
    '''
    
    ''' output: options = [
        {
            'group': 'spatial',
            'items': [
                {'label': 'E7.5.h5ad', 'value': '/data1/share/omics-viewer/stintev/datasets/public/spatial/E7.5.h5ad'}, 
                {'label': E7.75.h5ad', 'value': '/data1/share/omics-viewer/stintev/datasets/public/spatial/E7.75.h5ad'},
                {'label': E8.0.h5ad', 'value': '/data1/share/omics-viewer/stintev/datasets/public/spatial/E8.0.h5ad'},
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
                    dataset_dir = os.path.join(path_server_folder,'datasets',i['group'],current_user.id,dataset)
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
                                    if adata.endswith('h5ad') and os.path.exists( os.path.join(dataset_dir, adata) )
                                ]
                            }
                        )

    all_options = {
        'overview_panels': [ options for i in range(len(ctx.outputs_grouping['overview_panels'])) ],
        'LR_panels': options
    }
    
    return all_options

#endregion

#region global colormap

@dashapp.callback( # update global colormap
    Output('STORE_global_cmap-overview', 'data'),
    Input({'type': 'PlotPanel_item_select_column', 'index': ALL}, 'value'),
    State({'type': 'PlotPanel_item_select_info', 'index': ALL}, 'value'),
    State({'type': 'PlotPanel_item_select_sample', 'index': ALL}, 'value'),
    State('STORE_global_cmap-overview', 'data'),
    State('STORE_plotPanelsCurUUID-overview', 'data'),
)
def update_global_cmap(
    column: List[str], info: List[str], 
    path_sample: List[str], global_cmap: Dict,
    uuid_list: List[str]
):
    
    # 同时只会有一个PlotPanel触发，如果有多个触发，说明是linkage触发，则只处理第一个
    if isinstance(ctx.triggered_id, list):
        t_id = ctx.triggered_id[0]['index']
    else:
        t_id = ctx.triggered_id['index']
    t_index = uuid_list.index(t_id)
    
    adata = anndata.read_h5ad(
        path_sample[t_index], backed='r'
    )
    
    if info[t_index] == 'metadata' and column[t_index]:
        categories = adata.obs[column[t_index]].unique().tolist()
        global_cmap = assign_colors(categories, global_cmap)

        return global_cmap

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
    Output({'type': 'PlotPanel_store_curCategories', 'index': MATCH}, 'data'),
    Output({'type': 'PlotPanel_item_categoriesLegend', 'index': MATCH}, 'children'), # categoriesLegend in PlotPanel
    
    Input({'type': 'PlotPanel_item_select_column', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_column', 'index': MATCH}, 'id'), # get content of MATCH by id['index']
    Input({'type': 'PlotPanel_item_select_embedding', 'index': MATCH}, 'value'),
    
    Input({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_info', 'index': MATCH}, 'value'),
    Input({'type': 'DataFilter_store_preserved_cells', 'index': MATCH}, 'data'),
    Input('STORE_global_cmap-overview', 'data'),
    prevent_initial_call=True
)
def update_PlotPanel_figure(
    column: str, 
    id: dict, 
    embedding: str, 
    path_sample: str, 
    info: str, 
    preserved_cells: list[str],
    global_cmap: Dict
):
    
    adata = anndata.read_h5ad(
        path_sample, backed='r'
    )
    
    index = id['index']
    
    if preserved_cells is None or len(preserved_cells) < 1:
        preserved_cells = adata.obs_names
    
    if info == 'feature' and column and embedding and path_sample and info:
        figure = plot_feature_embedding(adata, preserved_cells, column, embedding)
        traceNumber = 1
        curCategories = None
        categoriesLegend = []
    elif info == 'metadata' and column and embedding and path_sample and info:
        figure = plot_metadata_embedding(adata, preserved_cells, column, embedding, global_cmap['cmap'])
        tmp_curCategories = adata[preserved_cells].obs[column].unique()
        traceNumber = len(tmp_curCategories)
        if adata.obs[column].dtypes == 'category':
            curCategories = [cat for cat in adata.obs[column].cat.categories.to_list() if cat in tmp_curCategories]
        else:
            curCategories = list(tmp_curCategories)
        categoriesLegend = PlotPanel.categoriesLegend(
            curCategories,  # order of px-legend, px-traces are all defined by series.cat.categories
            index = index, cmap = global_cmap['cmap']
        )
    else:
        raise PreventUpdate
    
    return figure, traceNumber, curCategories, categoriesLegend

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
    Output({'type': 'PanelLinkages_icon_linkage', 'index': ALL}, 'color'),
    inputs={
        'all_inputs': {
            'list_types': Input({'type': 'PanelLinkages_select_type', 'index': ALL}, 'value'),
            'list_panels': Input({'type': 'PanelLinkages_select_linkage', 'index': ALL}, 'value'),
            'list_apply': Input({'type': 'PanelLinkages_switch_apply', 'index': ALL}, 'checked'),
            'list_plotPanels_curUUID': Input('STORE_plotPanelsCurUUID-overview', 'data'), # triggered when add/delete PlotPanel
        }
    },
)
def apply_linkage(all_inputs):
    
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
                'property': 'value', 'value': ['column', '3D view'], 
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

    list_linkage_type = [item['value'] for item in ctx.args_grouping.all_inputs.list_types]
    list_selected_panels = [item['value'] for item in ctx.args_grouping.all_inputs.list_panels]
    list_linkage_apply = [item['value'] for item in ctx.args_grouping.all_inputs.list_apply]
    panels_curUUID = all_inputs['list_plotPanels_curUUID']
    color_sequence = ['red', 'blue', 'green']
    i = -1
    marks = [[]] * len(panels_curUUID) # panels marks
    colors = [] # linkage colors
    for type, panels, apply in zip(list_linkage_type, list_selected_panels, list_linkage_apply): # 对每条linkage
        # 颜色迭代器
        i += 1
        i = i % 3
        if apply and type: # 如果应用linkage
            colors.append(color_sequence[i])
            for panel in panels: # 对于选中的每个panel（uid）
                index = panels_curUUID.index(panel) # 找到uid对应的现有panel的位置
                marks[index] = marks[index] + [
                    PanelLinkages.linkage_mark(color=color_sequence[i])
                ] # 不能改成marks[index].append(...)或者marks[index]+=...，否则从marks为[[],[], ...]时每次都会会给所有panels添加mark
        else: # 不应用的话则灰色，不添加marks
            colors.append('gray')
    return marks, colors

#clientside
@dashapp.callback( # update linked panels
    output = {
        'samples': Output({'type': 'PlotPanel_item_select_sample', 'index': ALL}, 'value'),
        'embeddings': Output({'type': 'PlotPanel_item_select_embedding', 'index': ALL}, 'value'),
        'features' : Output({'type': 'PlotPanel_item_select_column', 'index': ALL}, 'value'),
        'infos': Output({'type': 'PlotPanel_item_select_info', 'index': ALL}, 'value'),
        'graphs' : Output({'type': 'PlotPanel_item_graph', 'index': ALL}, 'figure'),
        'highlight_index': Output({'type': 'PlotPanel_store_selectedCategories_index', 'index': ALL}, 'data'),
    },
    inputs = {
        'samples' : Input({'type': 'PlotPanel_item_select_sample', 'index': ALL}, 'value'),
        'embeddings': Input({'type': 'PlotPanel_item_select_embedding', 'index': ALL}, 'value'),
        'features' : Input({'type': 'PlotPanel_item_select_column', 'index': ALL}, 'value'),
        'infos': Input({'type': 'PlotPanel_item_select_info', 'index': ALL}, 'value'),
        'relayoutDatas' : Input({'type': 'PlotPanel_item_graph', 'index': ALL}, 'relayoutData'),
        'linkages_type': Input({'type': 'PanelLinkages_select_type', 'index': ALL}, 'value'),
        'linkages_panels': Input({'type': 'PanelLinkages_select_linkage', 'index': ALL}, 'value'),
        'linkages_apply': Input({'type': 'PanelLinkages_switch_apply', 'index': ALL}, 'checked'),
        'plotPanel_uuids': State('STORE_plotPanelsCurUUID-overview', 'data'),
        'curHighlight_index': Input({'type': 'PlotPanel_store_selectedCategories_index', 'index': ALL}, 'data'),
        'plotPanel_store_curCategories': State({'type': 'PlotPanel_store_curCategories', 'index': ALL}, 'data'),
    },
)
def update_linked_panels_sample(
    samples: List[str], embeddings: List[str], features: List[str], infos: List[str],
    relayoutDatas: List, 
    linkages_type: List[PanelLinkages.LinkageTypes],
    linkages_panels: List[str], # uuid
    linkages_apply: List[bool],
    plotPanel_uuids: List[str],
    plotPanel_store_curCategories: List[List[str]], # categories of each panel
    curHighlight_index: List[List[int]]
):
    
    return_dict = {}
    tid = ctx.triggered_id
    
    if tid and 'type' in tid and tid['type'].startswith('PanelLinkages_'): 
        # linkages热改动触发
        raise PreventUpdate # 暂时
    elif tid and 'type' in tid and tid['type'].startswith('PlotPanel_'): 
        # PlotPanel select改动触发
        if not any(linkages_apply): # 至少有linkage再执行
            raise PreventUpdate
        
        return_samples = [no_update]*len(plotPanel_uuids)
        return_embeddings = [no_update]*len(plotPanel_uuids)
        return_features = [no_update]*len(plotPanel_uuids)
        return_infos = [no_update]*len(plotPanel_uuids)
        return_graphs = [Patch()]*len(plotPanel_uuids)
        return_highlightIndex = [no_update]*len(plotPanel_uuids)
        
        panel_order = plotPanel_uuids.index(tid['index']) # 触发panel在队列中的位置
        
        if tid['type'] == 'PlotPanel_item_select_sample':
            for type, panels, apply in zip(linkages_type, linkages_panels, linkages_apply): # 对每条linkage
                if (apply is False) or (type is [None]) or (type is None): # 如果没有apply或者type没有选，跳过该条linkage
                    continue
                if tid['index'] not in panels: # 触发的panel没有在该条linkage被选中
                    continue
                
                if 'sample' in type and len(panels) >= 2: # 如果选项中有sample，且linkage选中的panel超过两个
                    sample_to_set = samples[panel_order]
                    for i,uid in enumerate(plotPanel_uuids):
                        if uid in panels and uid != tid['index']:
                            return_samples[i] = sample_to_set
            
        if tid['type'] == 'PlotPanel_item_select_embedding':  
            for type, panels, apply in zip(linkages_type, linkages_panels, linkages_apply): # 对每条linkage
                if (apply is False) or (type is [None]) or (type is None): # 如果没有apply或者type没有选，跳过该条linkage
                    continue
                if tid['index'] not in panels: # 触发的panel没有在该条linkage被选中
                    continue
                
                if 'embedding' in type and len(panels) >= 2: # embedding
                    embedding_to_set = embeddings[panel_order]
                    for i,uid in enumerate(plotPanel_uuids):
                        if uid in panels and uid != tid['index']:
                            return_embeddings[i] = embedding_to_set
                            
        if tid['type'] in ['PlotPanel_item_select_column','PlotPanel_item_select_info']:       
            for type, panels, apply in zip(linkages_type, linkages_panels, linkages_apply): # 对每条linkage
                if (apply is False) or (type is [None]) or (type is None): # 如果没有apply或者type没有选，跳过该条linkage
                    continue
                if tid['index'] not in panels: # 触发的panel没有在该条linkage被选中
                    continue
                    
                if 'feature' in type and len(panels) >= 2: # column
                    column_to_set = features[panel_order]
                    info_to_set = infos[panel_order]
                    for i,uid in enumerate(plotPanel_uuids):
                        if uid in panels and uid != tid['index']:
                            return_features[i] = column_to_set
                            return_infos[i] = info_to_set
                       
        if tid['type'] == 'PlotPanel_item_graph': 
            for type, panels, apply in zip(linkages_type, linkages_panels, linkages_apply): # 对每条linkage
                if (apply is False) or (type is [None]) or (type is None): # 如果没有apply或者type没有选，跳过该条linkage
                    continue
                if tid['index'] not in panels: # 触发的panel没有在该条linkage被选中
                    continue
                
                if '3D view' in type and len(panels) >= 2: # view
                    view_to_set = relayoutDatas[panel_order]
                    for i,uid in enumerate(plotPanel_uuids):
                        if uid in panels and uid != tid['index']:
                            if 'scene.camera' in view_to_set:
                                return_graphs[i]['layout']['scene']['camera'] = view_to_set['scene.camera']
                            if 'scene.aspectratio' in view_to_set:
                                return_graphs[i]['layout']['scene']['aspectmode'] = 'manual'
                                return_graphs[i]['layout']['scene']['aspectratio'] = view_to_set['scene.aspectratio']

        # categoriesLegend 触发
        if tid['type'] == 'PlotPanel_store_selectedCategories_index':
            for type, panels, apply in zip(linkages_type, linkages_panels, linkages_apply): # 对每条linkage
                if (apply is False) or (type is [None]) or (type is None): # 如果没有apply或者type没有选，跳过该条linkage
                    continue
                if tid['index'] not in panels: # 触发的panel没有在该条linkage被选中
                    continue
                
                if 'highlighting' in type and len(panels) >= 2: # highlighting
                    cat_A = plotPanel_store_curCategories[panel_order]
                    hlt_A = [ cat_A[i] for i in curHighlight_index[panel_order] ]
                    for i,uid in enumerate(plotPanel_uuids):
                        if uid in panels and uid != tid['index']:
                            cat_B = plotPanel_store_curCategories[i]
                            hlt_B = [ cat_B[i] for i in curHighlight_index[i] ]
                            hlt_after = set(hlt_A + hlt_B) - (set(cat_A)-set(hlt_A))
                            index_to_hlt = [ cat_B.index(cat) for cat in hlt_after if cat in cat_B]
                            return_highlightIndex[i] = index_to_hlt
            
            
            
        return_dict = {
            'samples': return_samples,
            'embeddings': return_embeddings,
            'features': return_features,
            'infos': return_infos,
            'graphs': return_graphs,
            'highlight_index': return_highlightIndex
        }
        
        return return_dict
    
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
    prevent_initial_call=True
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

#region ligand-receptor

#endregion
