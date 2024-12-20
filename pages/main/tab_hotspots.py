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

class TabHotspots:
    # const
    _width_sider = '15vw'
    _width_sider_collapsed = 60
    _width_drawer = 1000
    _rowHeight_plot_panel = 350
    _height_plot_panel_item: int = 300
    
    def __init__(self) -> None:
        
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
            label = 'Hotspots(CCC)',
            tab_id = 'TAB-Hotspots',
            children=[
                fac.Layout(
                    [
                        self.sider,
                        self.content,
                    ],
                    className='fac-Layout'
                ),
            ]
        )