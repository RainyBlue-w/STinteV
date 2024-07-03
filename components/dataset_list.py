import shutil
from typing import Literal
import dash_mantine_components as dmc
from dash_extensions.enrich import html
from dash_extensions.enrich import Output, Input, State, MATCH, ALL, callback, clientside_callback, no_update
from dash.exceptions import PreventUpdate
from dash import ctx, Patch, dcc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac

import os
import re

from flask_login import current_user
from numpy import place

from stintev.utils._io import read_txt
from stintev.config import PathConfig

class DatasetList:
    
    @staticmethod
    def alert_guest():
        return dmc.Alert(
            'You are browsing as a guest, sign in first to use this!',
            title = 'Permission Denied',
            color = 'violet',
            withCloseButton=False,
            icon = DashIconify(icon='mdi:denied', width=24),
            className='dmc-Alert-guest-dataset'
        )
        
    @staticmethod
    def alert_empty():
        return dmc.Alert(
            'There is no private dataset preserved, upload one first!',
            title = 'No dataset',
            color = 'blue',
            withCloseButton=False,
            icon = DashIconify(icon='iconoir:db-error', width=48),
            className='dmc-Alert-empty-dataset'
        )
    
    def grid_cols(self, group: Literal['public', 'private'], path_data_folder: str):
        content = []
        
        for dataset, files in self.datasets.items():
            content.append(
                dmc.GridCol(
                    id={'type': 'DatasetList-gridCol', 'index': f'{group}-{dataset}'},
                    span=3,
                    children = dmc.Card(
                        style={'height': '100%'},
                        radius='md',
                        shadow='sm',
                        withBorder=True,
                        children=[
                            dmc.Checkbox(
                                value=dataset,
                                size='lg',
                                label=[
                                    dmc.Stack(
                                        gap=0,
                                        children=[
                                            dmc.Group(
                                                children=[
                                                    dmc.Text(dataset, size='lg'),
                                                    dmc.Burger(
                                                        id={'type': 'DatasetList-burger-button', 'index': f'{group}-{dataset}'},
                                                        opened=False,
                                                        style={'position': 'absolute', 'right': '45px'}
                                                    ),
                                                    html.Div( # 避免popconfirm直接施加到position: absolute元素导致的悬浮层定位问题，absolute改为施加到该div上
                                                        children = fac.Popconfirm(
                                                            children = dmc.ActionIcon(
                                                                DashIconify(icon='gg:close-r', width=24),
                                                                size='lg', variant='light', color='red',
                                                            ),
                                                            id = {'type': 'DatasetList-close-confirm', 'index': f'{group}-{dataset}'},
                                                            title = 'Delete this dataset?',
                                                            className='fac-Popconfirm-delete-dataset',
                                                            locale='en-us',
                                                            trigger='click',
                                                            placement='bottom',
                                                            okText='Yes',
                                                            okButtonProps = dict(
                                                                size='middle', danger=True
                                                            ),
                                                            cancelText='No',
                                                            cancelButtonProps = dict(
                                                                size='middle', danger=False, type='ghost'
                                                            ),
                                                        ),
                                                        style={'position': 'absolute', 'right': '10px'}
                                                    )
                                                ] if group=='private' else [
                                                    dmc.Text(dataset, size='lg'),
                                                    dmc.Burger(
                                                        id={'type': 'DatasetList-burger-button', 'index': f'{group}-{dataset}'},
                                                        opened=False,
                                                        style={'position': 'absolute', 'right': '10px'}
                                                    )
                                                ],
                                            ),
                                            dmc.Space(h=5),
                                            html.Div(
                                                hidden=True,
                                                children=dmc.Stack(
                                                    [dmc.Text(file, size='md', c='gray') for file in files],
                                                    gap=0,
                                                ),
                                                id={'type': 'DatasetList-dataset-files', 'index': f'{group}-{dataset}'}
                                            ),
                                        ],
                                    )
                                ]
                            ),
                            dmc.Space(h=5),
                            dmc.CardSection(
                                children=[
                                    dmc.Text(
                                        children=read_txt(os.path.join(path_data_folder, dataset, 'description.txt')),
                                        c='gray', style={'margin-left': '10px'}
                                    )
                                ]
                            ),
                            dmc.Space(h=2),
                        ],
                    )
                )
            )
        
        if len(content) == 0:
            content.append(
                dmc.GridCol(
                    DatasetList.alert_empty(),
                    span=12
                )
            )
            
        return content
        
    def __init__(
        self, 
        path_data_folder: str, 
        group: str, # Public | Private
    ):
        
        if os.path.exists(path_data_folder) is False:
            os.makedirs(path_data_folder)
        
        self.datasets = {
            dataset : [
                file for file in
                sorted(os.listdir(os.path.join(path_data_folder, dataset)))
                if file.endswith('.h5ad')
            ]
            for dataset in 
            sorted(os.listdir(os.path.join(path_data_folder)))
        }
        
        self.list = dmc.CheckboxGroup(
            id={'type': 'DatasetList-checkboxGroup', 'index': group},
            value=[],
            children=[
                dmc.Grid(
                    id = {'type': 'DatasetList-checkboxGroup-Grid', 'index': group},
                    children=self.grid_cols(group=group, path_data_folder=path_data_folder),
                ),
                dcc.Store(
                    id={'type': 'DatasetList-Store-curIndexes', 'index': group},
                    data=[
                        f'{group}-{dataset}'
                        for dataset in self.datasets.keys()
                    ]
                ),
            ]
        )
        
clientside_callback( # 展开数据集文件列表
    """
    function(opened) {
        return !opened
    }
    """,
    Output({'type': 'DatasetList-dataset-files', 'index': MATCH}, 'hidden'),
    Input({'type': 'DatasetList-burger-button', 'index': MATCH}, 'opened')
)

@callback( # 删除private dataset
    Output('TABS_panel_private-dataset', 'children'),
    Output({'type': 'DatasetList-checkboxGroup-Grid', 'index': 'private'}, 'children'),
    Output({'type': 'DatasetList-Store-curIndexes', 'index': 'private'}, 'data'),
    
    Input({'type': 'DatasetList-close-confirm', 'index': ALL}, 'confirmCounts'),
    
    State({'type': 'DatasetList-Store-curIndexes', 'index': 'private'}, 'data'),
    prevent_initial_call=True
)
def delete_private_dataset(list_clicks, curIndexes):
    
    if all(list_clicks): # 有确认删除
    
        tid = ctx.triggered_id
        '''
        tid: {'index': 'private-test', 'type': 'DatasetList-close-button'}
        '''
        
        children_Grid = Patch()
        children_curIndexes = Patch()
        for i, index in enumerate(curIndexes):
            if index == tid['index']:
                del children_Grid[i]
                del children_curIndexes[i]
                dataset = re.findall(r'^\w+-(.+)',tid['index'])[0]
                shutil.rmtree(os.path.join(PathConfig.DATA_PATH, f'datasets/private/{current_user.username}/{dataset}'))
                
        if len(curIndexes) <= 1: # 删除最后一个dataset
            return DatasetList.alert_empty(), children_Grid, children_curIndexes
        else: # 删除后还有dataset
            return no_update, children_Grid, children_curIndexes
        
    else:   # 拦截生成组件时产生的触发
        raise PreventUpdate
            