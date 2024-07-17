import uuid
from dash_extensions.enrich import Output, Input, State, html, callback, clientside_callback, ClientsideFunction
from dash import dcc, ALL, MATCH, Patch, no_update, ctx
from dash.exceptions import PreventUpdate
from dash_iconify import DashIconify
import dash_mantine_components as dmc
import feffery_utils_components as fuc
import feffery_antd_components.alias as fac
import plotly.express as px
import plotly.graph_objects as go

from typing import Dict, Tuple, List
import uuid

from stintev.utils._io import *
from stintev.utils._plot import *
from .data_filter import DataFilter

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
        init_samples: List = None,
    ) -> None:
        '''
        index: str
            uuid to make unique
        '''
        self._index = index

        self._store = html.Div([
            dcc.Store(id = {'type': 'PlotPanel_store_traceNumber', 'index': self._index}, data=1),
            # store metadata columns for filtering
            dcc.Store(id = {'type': 'PlotPanel_store_metadata_columns', 'index': self._index}) 
        ])

        self._display_idx = dmc.Group(
            gap=3,
            children = [
                dmc.Group( # linkage marks
                    children = [],
                    id={'type': 'PlotPanel_linkage_marks', 'index': self._index},
                    gap=0,
                    style = {'position': 'absolute', 'left': '5px', 'top': '5px'}
                ),
                dmc.Badge( # panel idx badge
                    id={'type': 'PlotPanel_badge_panel_idx', 'index': self._index},
                    variant='light', color='gray',
                    style = {'position': 'absolute', 'right': '25px', 'top': '5px'}
                ),
            ],
        )

        self._data_button = fac.Popover(
            trigger = 'click',
            mouseLeaveDelay=0.1,
            placement = 'bottom',
            children = dmc.ActionIcon(
                children=DashIconify(icon='fluent:database-multiple-20-regular', width=24),
                size='md', variant='light', color='blue'
            ),
            content = dmc.Stack(
                [
                    dmc.Text('Sample', className='dmc-Text-select-label-PlotPanel'),
                    fac.Select(
                        placeholder='Choose a sample',
                        locale = 'en-us',
                        allowClear=False,
                        id = {'type': 'PlotPanel_item_select_sample', 'index': self._index},
                        options = init_samples if init_samples else []
                    ),
                    dmc.Space(h=3),
                    dmc.Text('Embedding', className='dmc-Text-select-label-PlotPanel'),
                    fac.Select(
                        placeholder='Choose an embedding',
                        locale = 'en-us',
                        allowClear=False,
                        id = {'type': 'PlotPanel_item_select_embedding', 'index': self._index}
                    ),
                    dmc.Space(h=3),
                    dmc.Text('Info', className='dmc-Text-select-label-PlotPanel'),
                    fac.Select(
                        placeholder='Info to plot',
                        locale = 'en-us',
                        allowClear=False,
                        id = {'type': 'PlotPanel_item_select_info', 'index': self._index},
                        options = [
                            {'label': 'feature', 'value': 'feature'},
                            {'label': 'metadata', 'value': 'metadata'},
                        ],
                    ),
                ],
                gap=0,
            )
        )

        self._settings_button = fac.Popover(
            trigger = 'click',
            mouseLeaveDelay=0.1,
            placement = 'bottom',
            children = dmc.ActionIcon(
                children=DashIconify(icon='fluent:settings-24-regular', width=24),
                size='md', variant='light', color='blue'
            ),
            content = dmc.Stack(
                [
                    dmc.Text('Point size', className='dmc-Text-select-label-PlotPanel'),
                    dmc.NumberInput(
                        id={'type': 'PlotPanel_item_pointSize', 'index': self._index},
                        value=2, step=0.5, min=0.5,
                    ),
                ],
                gap=0,
            )
        )

        self._filter_button = fac.Popover(
            trigger = 'click',
            mouseLeaveDelay=0.1,
            placement = 'bottom',
            children = dmc.ActionIcon(
                children=DashIconify(icon='iconoir:filter', width=24),
                size='md', variant='light', color='blue'
            ),
            content = dmc.Stack(
                gap=0,
                children=[
                    DataFilter(index=self._index).filter
                ]
            ),
            className='fac-Popover-filter'
        )

        self.grid_item_control = html.Div(
            dmc.Grid(
                columns=100,
                align='end',
                gutter='5px',
                children=[
                    dmc.GridCol(
                        children= dmc.Stack(
                            gap=3,
                            children = [
                                self._display_idx,
                                dmc.Group(
                                    gap=3,
                                    children=[
                                        self._data_button, 
                                        self._settings_button, 
                                        self._filter_button, 
                                        # self._linkage_button
                                    ]
                                ),
                            ]
                        ),
                        span='content'
                    ),
                    
                    # select column
                    dmc.GridCol(
                        dmc.Stack(
                            [
                                dmc.Text('Column', className='dmc-Text-select-label-PlotPanel'),
                                fac.Select(
                                    locale = 'en-us',
                                    allowClear=False,
                                    optionFilterMode = 'case-sensitive',
                                    id = {'type': 'PlotPanel_item_select_column', 'index': self._index},
                                ),
                            ],
                            gap=0
                        ),
                        span='auto'
                    ),
                    
                    # close button
                    dmc.GridCol(
                        dmc.ActionIcon(
                            DashIconify(icon='gg:close-r', width=24),
                            id = {'type': 'PlotPanel_item_button_delete', 'index': self._index},
                            size='md', variant='light', color='red'
                        ),
                        span='content'
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
                    # shadow = 'hover-shadow',
                    className = 'fuc-div-plotPanel-gridItem',
                    children=[
                        self.grid_item_control,
                        self._store,
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
                                'toImageButtonOptions': {'format': 'jpeg', 'scale': 2}
                            },
                            style = {'height': f'{self._height_plot_panel_item}px'}
                        )   
                    ]
                )
            ]
        )
