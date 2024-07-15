import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash import dcc
from dash_iconify import DashIconify

import os
from flask_login import current_user

from stintev.components import UploadPanel, DatasetList

class TabDataset:
    
    _width_sider = 400
    
    def __init__(self,  path_data_folder: str):

        # store for choosen datasets
        self.store_server_folder = dcc.Store(id = 'STORE_server_folder-dataset', data = path_data_folder)
        self.store_dataset = dcc.Store(id = 'STORE_choosen_dataset-dataset')
        
        self.tab = dbc.Tab(
            label = 'Dataset',
            tab_id = 'TAB-dataset',
            children=[
                self.store_server_folder,
                self.store_dataset,
                dmc.Tabs(
                    id='TABS_group-dataset',
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
                            children = [
                                DatasetList(
                                    path_data_folder=os.path.join(path_data_folder,'datasets','public'),
                                    group='public'
                                ).list,
                            ],
                            value = 'public'
                        ),
                        dmc.TabsPanel(
                            id = 'TABS_panel_private-dataset',
                            children = [
                                DatasetList(
                                    path_data_folder=os.path.join(path_data_folder,'datasets','private', current_user.username),
                                    group='private'
                                ).list if current_user.is_authenticated else 
                                DatasetList.alert_guest()
                            ],
                            value = 'private'
                        ),
                        dmc.TabsPanel(
                            id = 'TABS_panel_upload-dataset',
                            children = [
                                UploadPanel().panel
                            ],
                            value = 'upload',
                            className = 'dmc-TabsPanel-upload-dataset'
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
