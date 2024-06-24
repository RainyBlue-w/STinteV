from typing import Literal
from dash.dash_table.Format import Format, Group, Scheme, Symbol
import dash_extensions
from dash_extensions.enrich import html
import dash_extensions.enrich
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac
import feffery_utils_components as fuc
from dash_extensions.enrich import Output, Input, State, Serverside
from dash_extensions.enrich import callback, clientside_callback, ClientsideFunction, callback_context
from dash import dcc
from dash import ALL, MATCH, Patch, ctx
from dash.exceptions import PreventUpdate

import uuid
import os
import anndata

from stintev.page_templates._io.read import read_dataset
from stintev.page_templates.template_2d.dataset_list import DatasetList

class TabDataset:
    
    _width_sider = 400
    
    def __init__(self, path_server_folder: str):

        # store for choosen datasets
        self.store_server_folder = dcc.Store(id = 'STORE_server_folder-dataset', data = path_server_folder)
        self.store_dataset = dcc.Store(id = 'STORE_choosen_dataset-dataset')
        
        self.tab = dbc.Tab(
            label = 'Dataset',
            tab_id = 'TAB-dataset',
            children=[
                self.store_server_folder,
                self.store_dataset,
                dmc.Tabs(
                    orientation='vertical',
                    className='dmc-Tabs-dataset',
                    value='public',
                    children=[
                        dmc.TabsList(
                            [
                                self._tabs_tab(
                                    value='public',title="Public", 
                                    tips='Datasets hosted by the server'
                                ),
                                self._tabs_tab(
                                    value='private',title="Private",
                                    tips='Datasets uploaded by yourself'
                                ),
                                self._tabs_tab(
                                    value='upload',title="Upload",
                                    tips='Upload your own datasets'
                                ),
                            ],
                            className='dmc-TabsList-sider-dataset'
                        ),
                        dmc.TabsPanel(
                            id = 'TABS_panel_public-dataset',
                            children = [],
                            value = 'public'
                        ),
                        dmc.TabsPanel(
                            id = 'TABS_panel_private-dataset',
                            children = [],
                            value = 'private'
                        ),
                        dmc.TabsPanel(
                            id = 'TABS_panel_upload-dataset',
                            children = [],
                            value = 'upload'
                        ),
                    ],
                )
            ]
        )
        
    def _tabs_tab(self, value, title, tips):
        
        content = dmc.TabsTab(
            value=value,
            children=[
                dmc.Stack(
                    [
                        dmc.Text(title,size='lg'),
                        dmc.Text(tips,size='sm',c='gray'),
                    ],
                    gap='0px',
                    className='dmc-Stack-sider-dataset',
                )
            ],
            className='dmc-TabsTab-sider-dataset'
        )
        
        return content
    