from dash_extensions.enrich import Output, Input, State, Serverside
from dash_extensions.enrich import callback, clientside_callback, ClientsideFunction, callback_context
from dash import ALL, MATCH, Patch, ctx
from dash.exceptions import PreventUpdate
from .plot_panel import PlotPanel

import anndata
import os

from stintev.page_templates._plot import plot_feature_embedding, plot_metadata_embedding

    
clientside_callback( # open the drawer for setting panels
    """
    function(n_clicks){
        return true
    }
    """,
    Output('DRAWER_setting_panels-overview', 'opened'),
    Input('BUTTON_setting_panels-overview', 'n_clicks'),
)

@callback( # update the rowHeight for PlotPanel-items
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

@callback( # add & delete PlotPanel
   Output('FUCGRID_content-overview', 'children', allow_duplicate=True),
   Output('DRAWER_setting_panels_div-overview', 'children', allow_duplicate=True),
   Output('STORE_plotPanelsCurUUID-overview', 'data', allow_duplicate=True),
   Output('FUCGRID_content-overview', 'layouts'),
   Output('BUTTON_setting_panels_add-overview', 'disabled'),
   
   Input('BUTTON_setting_panels_add-overview', 'n_clicks'),
   Input({'type': 'PlotPanel_settings_button_delete', 'index': ALL}, 'n_clicks'),
   State('STORE_plotPanelsCurUUID-overview', 'data'),
   suppress_callback_exceptions=True,
   prevent_initial_call=True
)
def add_plot_panel(add, delete, uuid_list):
    
    import uuid
    
    tid = ctx.triggered_id
    children_grid = Patch()
    children_drawer = Patch()
    
    if tid == 'BUTTON_setting_panels_add-overview':
        if add and (len(uuid_list) <= 5):
            next_PlotPanel = PlotPanel(uuid.uuid1().hex)
            children_grid.append(
                next_PlotPanel.grid_item
            )
            children_drawer.append(
                next_PlotPanel.settings
            )
            uuid_list.append(
                next_PlotPanel._index
            )
            
    elif ('type' in tid) and (tid['type'] == 'PlotPanel_settings_button_delete') and (len(uuid_list) > 1):

        del_index = tid['index']
        for i, index in enumerate(uuid_list):
            if index == del_index:
                del children_grid[i]
                del children_drawer[i]
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
        
    return  children_grid, children_drawer, uuid_list, layouts, ban_button

# update feature options by search_value, rule: startswith
# @callback(
#     Output({'type': 'PlotPanel_item_select_column', 'index': MATCH}, 'options'),
#     Input({'type': 'PlotPanel_item_select_column', 'index': MATCH}, 'search_value'),
#     Input({'type': 'PlotPanel_item_select_type', 'index': MATCH}, 'value'),
# )
# def update_feature_options_startswith_search_value(search_value, type):
#     pass

#region update PlotPanel

@callback( # update dataset for tab_overview
    Output('STORE_choosen_dataset-dataset', 'data'),
    Input({'type': 'DatasetList-checkboxGroup', 'index': ALL}, 'value'),    
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

@callback( # load choosen datasets in tab_dataset

    Output({'type': 'PlotPanel_item_select_sample', 'index': ALL}, 'data'),

    Input('BUTTON_load_choosen_dataset-dataset', 'n_clicks'),
    State('STORE_choosen_dataset-dataset', 'data'),
    State('STORE_server_folder-dataset', 'data')
)
def load_choosen_datasets(n_clicks, choosen_dataset, path_server_folder):
    
    ''' choosen_dataset: [
        {'group': 'public', 'choosen': ['spatial']},
        {'group': 'private', 'choosen': ['test2']}
    ]
    '''
    
    ''' ctx.outputs_list: 
    [
        {'id': 'STORE_choosen_dataset_loaded-dataset', 'property': 'data'}, 
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
    if n_clicks:
        
        for i in choosen_dataset:
            if len(i['choosen']) > 0 : 
                for dataset in i['choosen']:
                    options.append(
                        {
                            'group': f'{dataset}-{i["group"]}',
                            'items': [
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
                    
        all_options = [ options for i in range(len(ctx.outputs_list[1]))]
        
        return all_options
    
    else:
        raise PreventUpdate


@callback( # update PlotPanel figure
    Output({'type': 'PlotPanel_item_graph', 'index': MATCH}, 'figure'),
    Input({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_type', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_column', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_embedding', 'index': MATCH}, 'value'),
)
def update_PlotPnel_figure(path_sample, type, column, embedding):
    
    adata = anndata.read_h5ad(
        path_sample, backed='r'
    ) 

    if type == 'feature':
        figure = plot_feature_embedding(adata, column, embedding)
    elif type == 'metadata':
        figure = plot_metadata_embedding(adata, column, embedding)
    else:
        raise PreventUpdate

    return figure

@callback( # update PlotPanel options for select_column
    Output({'type': 'PlotPanel_item_select_column', 'index': MATCH}, 'options'),
    Input({'type': 'PlotPanel_item_select_type', 'index': MATCH}, 'value'),
    Input({'type': 'PlotPanel_item_select_sample', 'index': MATCH}, 'value'),
)
def update_PlotPael_column_options(type, path_sample):

    adata = anndata.read_h5ad(
        path_sample, backed='r'
    )
    
    if type == 'feature':
        return [{'label': column, 'value': column} for column in adata.var_names]
    elif type == 'metadata':
        return [{'label': column, 'value': column} for column in adata.obs.columns]
    else:
        raise PreventUpdate
    
#endregion

