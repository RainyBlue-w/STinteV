import anndata
from dash_extensions.enrich import Output, Input, html, callback, clientside_callback, ClientsideFunction
from dash import dcc, ALL, MATCH, Patch
from dash_iconify import DashIconify
import dash_mantine_components as dmc
import feffery_utils_components as fuc
import plotly.express as px
import plotly.graph_objects as go

from typing import Dict, Tuple

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
        dataset: Dict[str, anndata.AnnData],
        adata: str,
        init_figure: str,
        init_figure_params: Dict
    ) -> None:
        '''
        index: str
            uuid to make unique
        '''
        self._index = index
        
        if init_figure == 'feature':
            self._figure = plot_feature_embedding(dataset[adata], **init_figure_params)
            options_column = sorted(
                [{'value': column, 'label': column} for column in dataset[adata].var_names],
                key = lambda x: x['label']
            )
            value_column = init_figure_params['feature']
        elif init_figure == 'metadata':
            self._figure = plot_metadata_embedding(dataset[adata], **init_figure_params)
            options_column = sorted(
                [{'value': column, 'label': column} for column in dataset[adata].obs.columns],
                key = lambda x: x['label']
            )
            value_column = init_figure_params['column']

        self.grid_item_control = html.Div(
            dmc.Grid(
                columns=100,
                align='end',
                gutter='5px',
                children=[
                    dmc.GridCol(
                        dmc.Select(
                            label = 'Sample',
                            size = 'xs',
                            id = {'type': 'PlotPanel_item_select_sample', 'index': self._index},
                            data = list(dataset.keys()),
                            value = adata
                        ),
                        span=30
                    ),
                    dmc.GridCol(
                        dmc.Select(
                            label = 'Type',
                            size = 'xs',
                            id = {'type': 'PlotPanel_item_select_type', 'index': self._index},
                            data = [
                                {'value': 'feature', 'label': 'feature'},
                                {'value': 'metadata', 'label': 'metadata'}
                            ],
                            value = init_figure
                        ),
                        span=25
                    ),
                    dmc.GridCol(
                        dmc.Stack(
                            [
                                dmc.Text('Column', className='dcc-DropDown-PlotPanel-item-label-column'),
                                dcc.Dropdown(
                                    options = options_column,
                                    value = value_column,
                                    id = {'type': 'PlotPanel_item_select_column', 'index': self._index},
                                    className='dcc-DropDown-PlotPanel-item-select-column',
                                    optionHeight=25,
                                    clearable=False
                                ),
                            ],
                            gap=0
                        ),
                        span=30
                    ),
                    dmc.GridCol(
                        dmc.Popover(
                            [
                                dmc.PopoverTarget(dmc.Button('...')),
                                dmc.PopoverDropdown(
                                    dmc.Stack(
                                        [
                                            dmc.Select(
                                                label = 'Plot type',
                                                data = ['Embedding'],
                                                value = 'Embedding'
                                            ),
                                            dmc.Select(
                                                label = 'Embedding',
                                                id = {'type': 'PlotPanel_item_select_embedding', 'index': self._index},
                                                data = list(dataset[adata].obsm.keys()),
                                                value = init_figure_params['embedding'],
                                                selectFirstOptionOnChange=True
                                            )
                                        ]
                                    )
                                ),
                            ]
                        ),
                        span=15
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
                        dcc.Graph(
                            id={'type': 'PlotPanel_item_graph', 'index': self._index},
                            figure=self._figure,
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
