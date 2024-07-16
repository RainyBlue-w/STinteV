from dash_extensions.enrich import Output, Input, State, html, callback, clientside_callback, ClientsideFunction
from dash import dcc, ALL, MATCH, Patch, no_update, ctx
from dash.exceptions import PreventUpdate
from dash_iconify import DashIconify
import dash_mantine_components as dmc
import feffery_utils_components as fuc
import feffery_antd_components.alias as fac
import plotly.express as px
import plotly.graph_objects as go

from typing import Dict, Literal, Tuple, List

class DataFilter:
    
    def __init__(self, index):
        self._index = index
        
        self.store_preserved_cells = dcc.Store(id={'type': 'DataFilter_store_preserved_cells', 'index': self._index})
        
        self.filter = html.Div(
            className='div-DataFilter',
            children=[
                self.store_preserved_cells,
                dmc.Grid(
                    align='center',
                    children=[
                        dmc.GridCol(
                            [
                            dmc.Text('Column', className='dmc-Text-select-label-PlotPanel'),
                            dcc.Dropdown(id={'type': 'DataFilter_select_column', 'index': self._index}, clearable=True, searchable=True,
                                        persistence=True, persistence_type='local'),
                            ],
                            span=7
                        ),
                        dmc.GridCol(
                            [
                            dmc.Text('Type', className='dmc-Text-select-label-PlotPanel'),
                            dcc.Dropdown(['numeric', 'categorical'], id={'type': 'DataFilter_select_type', 'index': self._index}, 
                                        clearable=False, searchable=False, persistence=True, persistence_type='local'),
                            ],
                            span=5
                        ),
                        dmc.GridCol(
                            children=[],
                            id={'type': 'DataFilter_filter_body', 'index': self._index}, 
                            span=12,
                        ),
                        dmc.GridCol(
                            dmc.Switch(
                                id={'type': 'DataFilter_switch_apply', 'index': self._index}, 
                                onLabel='ON', offLabel='OFF', 
                                size='lg', checked=False
                            ), 
                            span=3
                        ),
                        dmc.GridCol(
                            dmc.Text('Selected cells: ', id={'type': 'DataFilter_text_number', 'index': self._index}, className='dmc-Text-sidebar-tips'),
                            span=9
                        )
                    ]
                ),
                # dmc.Space(h=5),
                # dmc.Divider(variant='dashed'),
                # dmc.Space(h=5),
            ]
        )
    
    @staticmethod
    def numeric_filter(index, column, min=None, max=None):
        return html.Div(
            className='DataFilter-div-filterBody-numeric',
            children=[
                dmc.Grid(
                    [
                        dmc.GridCol(dmc.Text(column), span='content'),
                        dmc.GridCol(dmc.Text(' ≥ '), span='content'),
                        dmc.GridCol(
                            dmc.NumberInput(value=min, decimalScale=4, step=0.1, id={'type': 'DataFilter_numberInput_left', 'index': index}),
                            span='content'
                        ),
                        dmc.GridCol(dmc.Text(column), span='content'),
                        dmc.GridCol(dmc.Text(' ≤ '), span='content'),
                        dmc.GridCol(
                            dmc.NumberInput(value=max, decimalScale=4, step=0.1, id={'type': 'DataFilter_numberInput_right', 'index': index}),
                            span='content'
                        ),
                    ],
                ),
                # solve nonexistent id problem(callback)
                html.Div(id={'type': 'DataFilter_transfer', 'index': index})
            ],
        )
    
    @staticmethod
    def categorical_filter(index, options):
        return html.Div(
            className='DataFilter-div-filterBody-categorical',
            children=[
                fac.Transfer(
                    dataSource=[ {'key': i, 'title': i} for i in options ],
                    locale='en-us',
                    targetKeys=[],
                    id={'type': 'DataFilter_transfer', 'index': index}
                ),
                # solve nonexistent id problem(callback)
                html.Div(id={'type': 'DataFilter_numberInput_left', 'index': index}),
                html.Div(id={'type': 'DataFilter_numberInput_right', 'index': index}),
            ]
        )