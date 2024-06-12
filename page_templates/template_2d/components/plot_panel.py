from curses import panel
from dash_extensions.enrich import Output, Input, html, callback, clientside_callback, ClientsideFunction
from dash import dcc, ALL, MATCH, Patch
from dash_iconify import DashIconify
import dash_mantine_components as dmc
import feffery_utils_components as fuc
import plotly.express as px
import plotly.graph_objects as go

class PlotPanel:
    
    '''
    class to generate a plot panel, including the graph and the settings
    '''
    
    # const
    _index = None
    _width_drawer_card = 220
    _height_plot_panel_item = 300
    
    # widgets
    grid_item=None
    settings=None
    
    def __init__(self, index: str) -> None:
        '''
        index: str
            uuid to make unique
        '''
        
        self._index = index
        
        self.grid_item = fuc.FefferyGridItem(
            key=str(self._index),
            id={'type': 'PLOTPANEL_item-2d', 'index': self._index},
            className = 'fuc-GirdItem',
            children=[
                fuc.FefferyDiv(
                    id = {'type':'PLOTPANEL_item_div-2d', 'index': self._index},
                    shadow = 'hover-shadow',
                    children=[
                        dcc.Graph(
                            id={'type': 'PLOTPANEL_item_graph-2d', 'index': self._index},
                            figure=px.scatter(x=[1], y=[1]).update_layout(
                                    # plot_bgcolor = '#ffffff', 
                                    uirevision='constant',
                                    coloraxis = {
                                        'colorbar' : {'tickformat': '4.2f'}
                                    },
                                    margin=dict(l=0, r=0, t=0, b=0), autosize=True
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
        
        self.settings = fuc.FefferyDiv(
            id = {'type':'PLOTPANEL_settings_card_div-2d', 'index': self._index},
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
                                            id = {'type': 'PLOTPANEL_settings_select_sample-2d', 'index': self._index},
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
                                            id = {'type': 'PLOTPANEL_settings_select_column-2d', 'index': self._index},
                                            data = []
                                        )    
                                    ]
                                ),
                                dmc.GridCol(
                                    span=100,
                                    children=[
                                        dmc.Button(
                                            'Delete panel', id={'type': 'PLOTPANEL_settings_button_delete-2d', 'index': self._index},
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
    Output({'type':'PLOTPANEL_settings_card_div-2d', 'index': MATCH}, 'shadow'),
    Output({'type':'PLOTPANEL_item_div-2d', 'index': MATCH}, 'shadow'),
    Input({'type':'PLOTPANEL_settings_card_div-2d', 'index': MATCH}, 'isHovering'),
    Input({'type':'PLOTPANEL_item_div-2d', 'index': MATCH}, 'isHovering'),
)

