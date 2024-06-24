import anndata
from dash_extensions.enrich import Output, Input, html, callback, clientside_callback, ClientsideFunction
from dash import dcc, ALL, MATCH, Patch
from dash_iconify import DashIconify
import dash_mantine_components as dmc
import feffery_utils_components as fuc
import feffery_antd_components.alias as fac
import plotly.express as px
import plotly.graph_objects as go

from typing import Dict, Tuple, List

from stintev.page_templates._io import *
from stintev.page_templates._plot import *

class PlotPanel:
    
    '''
    class to generate a plot panel, including the graph and the settings
    '''

    # const
    _index: str = None
    _width_drawer_card: int = 220
    _height_plot_panel_item: int = 300
    
    # widgets
    grid_item: fuc.FefferyGridItem = None
    settings: fuc.FefferyDiv = None

    # figure
    _figure: go.Figure = None
    
    def __init__(
        self, 
        index: str, 
        init_samples: List = None
    ) -> None:
        '''
        index: str
            uuid to make unique
        '''
        self._index = index

        self.grid_item_control = html.Div(
            dmc.Grid(
                columns=100,
                align='end',
                gutter='5px',
                children=[
                    dmc.GridCol(
                        children = fac.Popover(
                            mouseLeaveDelay=0.1,
                            placement = 'bottom',
                            children = dmc.ActionIcon(
                                children=DashIconify(icon='fluent:database-multiple-20-regular', width=24),
                                size='lg', variant='light', color='blue'
                            ),
                            content = dmc.Stack(
                                [
                                    dmc.Text('Sample', className='dcc-DropDown-PlotPanel-item-label-column'),
                                    fac.Select(
                                        placeholder='Choose a sample',
                                        locale = 'en-us',
                                        allowClear=False,
                                        id = {'type': 'PlotPanel_item_select_sample', 'index': self._index},
                                        options = init_samples if init_samples else []
                                    ),
                                    dmc.Space(h=3),
                                    dmc.Text('Embedding', className='dcc-DropDown-PlotPanel-item-label-column'),
                                    fac.Select(
                                        placeholder='Choose an embedding',
                                        locale = 'en-us',
                                        allowClear=False,
                                        id = {'type': 'PlotPanel_item_select_embedding', 'index': self._index}
                                    ),
                                ],
                                gap=0,
                            )
                        ),
                        span=9
                    ),
                    dmc.GridCol(
                       children = fac.Popover(
                            mouseLeaveDelay=0.1,
                            placement = 'bottom',
                            children = dmc.ActionIcon(
                                children=DashIconify(icon='fluent:settings-24-regular', width=24),
                                size='lg', variant='light', color='blue'
                            ),
                            content = dmc.Stack(
                                [
                                    dmc.Text('Point size', className='dcc-DropDown-PlotPanel-item-label-column'),
                                    dmc.NumberInput(
                                        id={'type': 'PlotPanel_item_pointSize', 'index': self._index},
                                        value=2, step=0.5, min=0.5,
                                    ),
                                ],
                                gap=0,
                            )
                        ),
                        span=9
                    ),
                    dmc.GridCol(
                        dmc.Stack(
                            [
                                dmc.Text('Info', className='dcc-DropDown-PlotPanel-item-label-column'),
                                fac.Select(
                                    locale = 'en-us',
                                    allowClear=False,
                                    id = {'type': 'PlotPanel_item_select_info', 'index': self._index},
                                    options = [
                                        {'label': 'feature', 'value': 'feature'},
                                        {'label': 'metadata', 'value': 'metadata'},
                                    ],
                                ),
                            ],
                            gap=0
                        ),
                        span=25
                    ),
                    dmc.GridCol(
                        dmc.Stack(
                            [
                                dmc.Text('Column', className='dcc-DropDown-PlotPanel-item-label-column'),
                                fac.Select(
                                    locale = 'en-us',
                                    allowClear=False,
                                    optionFilterMode = 'case-sensitive',
                                    id = {'type': 'PlotPanel_item_select_column', 'index': self._index}
                                ),
                            ],
                            gap=0
                        ),
                        span=40
                    ),
                    dmc.GridCol(
                        dmc.ActionIcon(
                            DashIconify(icon='gg:close-r', width=24),
                            size='lg', variant='light', color='red'
                        ),
                        span=10
                    )
                ]
            )
        )

        self.grid_item = fuc.FefferyGridItem(
            key=str(self._index),
            id={'type': 'PlotPanel_item', 'index': self._index},
            className = 'fuc-GirdItem',
            children=[
                fuc.FefferyDiv(
                    id = {'type':'PlotPanel_item_div', 'index': self._index},
                    shadow = 'hover-shadow',
                    className = 'fuc-div-plotPanel-gridItem',
                    children=[
                        self.grid_item_control,
                        dcc.Store(id = {'type': 'PlotPanel_store_traceNumber', 'index': self._index}, data=1),
                        dcc.Graph(
                            id={'type': 'PlotPanel_item_graph', 'index': self._index},
                            figure=px.scatter().update_layout(
                                margin=dict(l=0, r=0, t=0, b=0),
                                plot_bgcolor = '#ffffff', 
                                uirevision='constant',
                                legend_itemsizing = 'constant'
                            ).update_xaxes(visible=False).update_yaxes(visible=False),
                            responsive = True,
                            config = {
                                'autosizable': True,
                            },
                            style = {'height': f'{self._height_plot_panel_item}px'}
                        )   
                    ]
                )
            ]
        )
        
        self.sider_settings = fuc.FefferyDiv(
            id = {'type':'PlotPanel_settings_card_div', 'index': self._index},
            shadow = 'hover-shadow',
            className = 'fuc-Div-plotPanel-settings',
            children = [
                dmc.Card(
                    withBorder = True,
                    shadow = None,
                    radius = "md",
                    w = self._width_drawer_card,
                    children = [
                        dmc.Grid(
                            columns=100,
                            children=[
                                dmc.GridCol(
                                    span=52,
                                    children=[
                                        dmc.Select(
                                            label = 'Sample',
                                            # placeholder='sample',
                                            id = {'type': 'PlotPanel_settings_select_sample', 'index': self._index},
                                            data = [
                                                {'value': 'E7.5', 'label': 'E7.5'},
                                                {'value': 'E7.75', 'label': 'E7.75'},
                                                {'value': 'E8.0', 'label': 'E8.0'}
                                            ],
                                            value = 'E7.5'
                                        )
                                    ]
                                ),
                                dmc.GridCol(
                                    span=48,
                                    children=[
                                        dmc.Select(
                                            label = 'Type',
                                            # placeholder='Exp/Obs',
                                            data = [
                                                {'value': 'Exp', 'label': 'Exp'},
                                                {'value': 'Obs', 'label': 'Obs'}
                                            ],
                                            value = 'Exp'
                                        )    
                                    ]
                                ),
                                dmc.GridCol(
                                    span=100,
                                    children=[
                                        dmc.Select(
                                            label = 'Column',
                                            id = {'type': 'PlotPanel_settings_select_column', 'index': self._index},
                                            data = []
                                        )    
                                    ]
                                ),
                                dmc.GridCol(
                                    span=100,
                                    children=[
                                        dmc.Button(
                                            'Delete panel', id={'type': 'PlotPanel_settings_button_delete', 'index': self._index},
                                            color='red', fullWidth=True, variant="subtle",
                                            leftSection=DashIconify(icon="fluent:subtract-square-20-regular", width=20),
                                        )
                                    ]
                                )
                            ]
                        )
                    ]
                )
            ],
        )


clientside_callback(
    ClientsideFunction(
        namespace='plot_panel_2d',
        function_name='plotpanel_sync_shadow_between_card_and_item'
    ),
    Output({'type':'PlotPanel_settings_card_div', 'index': MATCH}, 'shadow'),
    Output({'type':'PlotPanel_item_div', 'index': MATCH}, 'shadow'),
    Input({'type':'PlotPanel_settings_card_div', 'index': MATCH}, 'isHovering'),
    Input({'type':'PlotPanel_item_div', 'index': MATCH}, 'isHovering'),
)
