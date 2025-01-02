
from dash_extensions.enrich import Output, Input, State, html, callback, clientside_callback, ClientsideFunction
from dash import dcc, ALL, MATCH, Patch, no_update, ctx, set_props
from dash.exceptions import PreventUpdate
from dash_iconify import DashIconify
import dash_mantine_components as dmc
import feffery_utils_components as fuc
import feffery_antd_components.alias as fac
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go

from typing import Dict, Tuple, List
import uuid
import re

from stintev.utils._io import *
from stintev.utils._plot import *
from .data_filter import DataFilter

class PlotPanel:

    '''
    class to generate a plot panel, including the graph and the settings
    '''

    # const
    _index: str | None = None
    _width_drawer_card: int = 220
    _height_plot_panel_item: int = 300

    # widgets
    grid_item: fuc.FefferyGridItem = None
    settings: fuc.FefferyDiv = None

    # figure
    _figure: go.Figure = None

    @staticmethod
    def categoriesLegend(
        categories: List[str], 
        index: str,
        cmap: Dict[str, str] | None, 
    ):
        if cmap is None:
            seq_colors = px.colors.qualitative.Alphabet
            target_len = len(categories)
            seq_colors = seq_colors*(target_len // len(seq_colors) + 1)
            seq_colors = seq_colors[0:target_len]
            cmap = dict(zip(categories, seq_colors))
            
        chips = dmc.Stack(
            [
                dmc.Grid(                    
                    gutter=2,
                    children=[
                        dmc.GridCol(
                            fac.Button(
                                icon=DashIconify(icon='system-uicons:reverse', width=16),
                                variant='filled', color='default', block=True,
                                id = {'type': 'PlotPanel_item_controlChip_invert', 'index': index},
                            ),
                            span=4
                        ),
                        dmc.GridCol(
                            fac.Button(
                                icon=DashIconify(icon='fluent:border-none-20-regular', width=16),
                                variant='filled', color='default', block=True,
                                id = {'type': 'PlotPanel_item_controlChip_clear', 'index': index},
                            ),
                            span=4,
                        ),
                        dmc.GridCol(
                            fac.Button(
                                icon=DashIconify(icon='fluent:checkbox-indeterminate-20-regular', width=16),
                                variant='filled', color='default', block=True,
                                id = {'type': 'PlotPanel_item_controlChip_all', 'index': index},
                            ),
                            span=4
                        ),
                        html.Div(
                            [
                                dbc.Tooltip( i.capitalize(), target={'type': f'PlotPanel_item_controlChip_{i}', 'index': index}, placement='top')
                                for i in ['invert', 'clear', 'all']
                            ],
                        ),
                    ]
                ),
                dmc.ChipGroup(
                    id = {'type': 'PlotPanel_item_categoriesChipGroup', 'index': index},
                    multiple=True,
                    deselectable=True,
                    value = categories,
                    children=[
                        dmc.Chip(
                            children=cat, value=cat, checked=True, color=cmap[cat],
                            size='xs', variant='filled', type='radio', autoContrast=True,
                            classNames = {'label': 'PlotPanel-legend-chip-label'},
                            id = {'type': 'PlotPanel_item_categoriesChip', 'index': f'{index}-{cat}'}
                        ) 
                        for cat in categories
                    ],
                ),
            ],
            justify = 'flex-start',
            gap=1,
        )
        return chips

    def __init__(
        self, 
        index: str, 
        init_samples: List | None = None,
    ) -> None:
        '''
        index: str
            uuid to make unique
        '''
        self._index = index

        self._store = html.Div([
            dcc.Store(id = {'type': 'PlotPanel_store_traceNumber', 'index': self._index}, data=1),
            # store current categories (in order) for trace-highlighting
            dcc.Store(id = {'type': 'PlotPanel_store_curCategories', 'index': self._index}, data=None), 
            # store selected categories (index)
            dcc.Store(id = {'type': 'PlotPanel_store_selectedCategories_index', 'index': self._index}, data=None),
            # store metadata columns for filtering
            dcc.Store(id = {'type': 'PlotPanel_store_metadata_columns', 'index': self._index}),
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
                                dmc.Text('Feature', className='dmc-Text-select-label-PlotPanel'),
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

        self.grid_item_graph = html.Div(
            dmc.Grid(
                gutter=3,
                children = [
                    # figure
                    dmc.GridCol(
                        dcc.Graph(
                            id={'type': 'PlotPanel_item_graph', 'index': self._index},
                            figure=px.scatter().update_layout(
                                margin=dict(l=0, r=0, t=0, b=0),
                                plot_bgcolor = '#ffffff', 
                                uirevision='constant',
                                legend_itemsizing = 'constant',
                            ).update_xaxes(visible=False).update_yaxes(visible=False),
                            responsive = True,
                            config = {
                                'autosizable': True,
                                'toImageButtonOptions': {'format': 'jpeg', 'scale': 2},
                                'showTips': True
                            },
                            style = {'height': f'{self._height_plot_panel_item}px'},
                        ),
                        id = {'type': 'PlotPanel_item_graph_leftCol', 'index': self._index},
                        span=9
                    ),
                    # legend
                    dmc.GridCol(
                        html.Div(
                            id={'type': 'PlotPanel_item_categoriesLegend', 'index': self._index}, 
                            style={'height': f'{self._height_plot_panel_item}px', 'overflow-y': 'auto'}
                        ),
                        id = {'type': 'PlotPanel_item_graph_rightCol', 'index': self._index},
                        span=3
                    )
                ],
            )
        )
        
        self.grid_item = fuc.FefferyGridItem(
            key=str(self._index),
            id={'type': 'PlotPanel_item', 'index': self._index},
            className = 'fuc-GirdItem',
            children=[
                fuc.FefferyDiv(
                    id = {'type':'PlotPanel_item_div', 'index': self._index},
                    className = 'fuc-div-plotPanel-gridItem',
                    children = dmc.Stack(
                        gap=0,
                        children=[
                            self.grid_item_control,
                            self._store,
                            self.grid_item_graph,
                        ]
                    )
                )
            ]
        )
        
#clientside
@callback( # update store_selectedCatergories_index on categoriesChipGroup changed
    Output({'type': 'PlotPanel_store_selectedCategories_index', 'index': MATCH}, 'data'),
    Output({'type': 'PlotPanel_item_graph', 'index': MATCH}, 'figure'),
    
    Input({'type': 'PlotPanel_item_categoriesChipGroup', 'index': MATCH}, 'value'),
    
    State({'type': 'PlotPanel_store_curCategories', 'index': MATCH}, 'data'),
    State({'type': 'PlotPanel_item_pointSize', 'index': MATCH}, 'value'), # cur point size

)
def update_figure_traces_highlighting_on_catLegend_clicking(categories_selected, curCategories, pt_size):

    return_index = [curCategories.index(cat) for cat in categories_selected]
    
    patch_fig = Patch()
    for i in range(len(curCategories)):
        if curCategories[i] in categories_selected:
            patch_fig['data'][i]['marker']['size'] = pt_size
        else:
            patch_fig['data'][i]['marker']['size'] = pt_size/5

    return return_index, patch_fig

@callback( # update categoriesChipGroup on store_selectedCategories_index changed 
    Output({'type': 'PlotPanel_item_graph', 'index': MATCH}, 'figure'),
    Output({'type': 'PlotPanel_item_categoriesChipGroup', 'index': MATCH}, 'value'),
    
    Input({'type': 'PlotPanel_store_selectedCategories_index', 'index': MATCH}, 'data'),
    State({'type': 'PlotPanel_store_curCategories', 'index': MATCH}, 'data'),
    State({'type': 'PlotPanel_item_pointSize', 'index': MATCH}, 'value'), # cur point size
)
def update_catLegend_on_store_selectedCategories_index_change(selected_index, curCategories, pt_size):
    patch_fig = Patch()
    for i in range(len(curCategories)):
        if i in selected_index:
            patch_fig['data'][i]['marker']['size'] = pt_size
        else:
            patch_fig['data'][i]['marker']['size'] = pt_size/5
    return patch_fig, [curCategories[idx] for idx in selected_index]

@callback(  # legend invert, clear, all
    Output({'type': 'PlotPanel_item_categoriesChipGroup', 'index': MATCH}, 'value'),

    Input({'type': 'PlotPanel_item_controlChip_invert', 'index': MATCH}, 'nClicks'),
    Input({'type': 'PlotPanel_item_controlChip_clear', 'index': MATCH}, 'nClicks'),
    Input({'type': 'PlotPanel_item_controlChip_all', 'index': MATCH}, 'nClicks'),
    
    State({'type': 'PlotPanel_item_categoriesChipGroup', 'index': MATCH}, 'value'),
    State({'type': 'PlotPanel_store_curCategories', 'index': MATCH}, 'data'),
)
def controlChip_invert_clear_all(nClicks_invert, nClicks_clear, nClicks_all, selected_categories, curCategories):
    
    tid = ctx.triggered_id
    
    if tid['type'] == 'PlotPanel_item_controlChip_invert' and nClicks_invert:
        return [cat for cat in curCategories if cat not in selected_categories]
    
    elif tid['type'] == 'PlotPanel_item_controlChip_clear' and nClicks_clear:
        return []
    
    elif tid['type'] == 'PlotPanel_item_controlChip_all' and nClicks_all:
        return curCategories
    
    raise PreventUpdate

# clientside
@callback( # update GridCol span on info type (feature/metadata) changed
    Output({'type': 'PlotPanel_item_graph_leftCol', 'index': MATCH}, 'span'),
    Output({'type': 'PlotPanel_item_graph_rightCol', 'index': MATCH}, 'span'),
    Input({'type': 'PlotPanel_item_select_info', 'index': MATCH}, 'value'),
)
def update_gridCol_span_on_infoType_changed(info_type):
    
    if info_type == 'feature':
        return 'auto', 'content'
    elif info_type == 'metadata':
        return 9, 3

    raise PreventUpdate