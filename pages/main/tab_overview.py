from os import name
from dash_extensions.enrich import html
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac
import feffery_utils_components as fuc
from dash import dcc
from dash_extensions.enrich import clientside_callback, Output, Input, ClientsideFunction

import uuid
from flask_login import current_user

from stintev.utils._plot import *
from stintev.components import PlotPanel, PanelLinkages

class TabOverview():
    
    # const
    _width_sider = '15vw'
    _width_sider_collapsed = 60
    _width_drawer = 250
    _rowHeight_plot_panel = 350
    
    # widgets
    store_dataset = None
    control_plot_settings = None
    
    control_plot_panels = None
    drawer_plot_panels = None
    _n_init_plot_panels = 2
    _init_plot_panels = None
    
    sider = None
    content = None
    tab = None
    
    # data
    dataset = None
    
    # Tab init
    def __init__(self, path_data_folder: str) -> None:

        self._init_plot_panels = [ 
            PlotPanel(index=uuid.uuid1().hex),
            PlotPanel(index=uuid.uuid1().hex),
        ]

        self._store = html.Div(
            [
                dcc.Store(
                    id='STORE_plotPanelsCurUUID-overview', 
                    data = [self._init_plot_panels[i]._index for i in range( self._n_init_plot_panels )]
                ),
                dcc.Store(
                    id = 'STORE_global_cmap-overview', storage_type = 'local',
                    data = {'i': 0, 'cmap': {}}
                ),
            ]
        )

        self.drawer_panel_linkage = None

        # accordion item in sider (Plot settings)
        self.control_plot_settings = dmc.AccordionItem(
            children=[
                dmc.AccordionControl(
                    dmc.Text('Plot settings', className='dmc-Text-accordionControl'),
                ),
                dmc.AccordionPanel(
                    children=[
                        # ----Points----
                        html.Div([
                            dmc.Divider(label = 'Points', labelPosition='center', variant='dashed', className='dmc-divider-sidebar-inline'),
                            dmc.Grid(
                                [
                                    dmc.GridCol(
                                        dmc.NumberInput(
                                            label = 'Point size',
                                            description = 'Size of the points in all the panels',
                                            id='NUMBERINPUT_scatter3dPointsize_3D',
                                            value=2, step=0.5, min=0.1,
                                            persistence='local',
                                        ),
                                        span=12
                                    ),
                                ],
                            ),
                            dmc.Space(h=5),        
                        ]),
                    ]
                )
            ],
            value='control_plot_settings'
        )

        # accordion item in sider (Plot panels)
        self.control_plot_panels = dmc.AccordionItem(
            children=[
                dmc.AccordionControl(
                    dmc.Text('Plot panels', className='dmc-Text-accordionControl'),
                ),
                dmc.AccordionPanel(
                    children = [
                        dmc.Grid(
                            children = [
                                dmc.GridCol(
                                    children = [
                                        dmc.NumberInput(
                                            label = 'Row height',
                                            id = 'NUMBERINPUT_setting_panels-overview',
                                            description = 'The minimal height of plot panel',
                                            value=self._rowHeight_plot_panel,
                                            min=200,
                                            step=25,
                                            persistence = True,
                                            persistence_type = 'local',
                                        )
                                    ]
                                ),
                                dmc.GridCol(
                                    children = [
                                        dmc.Button(
                                            'Add a panel', fullWidth=True, id='BUTTON_setting_panels_add-overview',
                                            leftSection=DashIconify(icon='fluent:add-square-20-regular', width=20),
                                        ),
                                    ],
                                    span=12
                                )
                            ]
                        ),
                    ]
                ),
            ],
            value='control_plot_panels'
        )

        # accordion item in sider (Linkage)
        self.control_linkage = dmc.AccordionItem(
            children=[
                dmc.AccordionControl(
                    dmc.Text('Panel linkages', className='dmc-Text-accordionControl'),
                ),
                dmc.AccordionPanel(
                    PanelLinkages.linkage_settings()
                )
            ],
            value = 'control_linkage'
        )

        self.control_colors = dmc.AccordionItem(
            children = [
                dmc.AccordionControl(
                    dmc.Text('Colors', className='dmc-Text-accordionControl'),
                ),
                dmc.AccordionPanel(
                    children=[
                        dmc.Button(
                            'Color settings', fullWidth=True, id='BUTTON_color_settings-overview',
                            leftSection=DashIconify(icon='proicons:color-palette', width=24),
                        )
                    ]
                ),
            ],
            value = 'Colors'
        )
        
        # sider in left
        self.sider = fac.Sider(
            collapsible = False,
            width = self._width_sider,
            theme = 'light',
            children = fac.Affix(
                html.Div([
                    dmc.Badge(
                        'Settings', color='blue', variant='light', 
                        radius='xs', size='xl', fullWidth=True,
                        leftSection=DashIconify(icon='fluent:settings-20-regular', width=20)
                    ),
                    dmc.Accordion(
                        multiple=True,
                        chevronPosition='right',
                        variant='seperated',
                        children=[
                            self.control_plot_settings,
                            self.control_plot_panels,
                            self.control_linkage,
                        ],
                        value=['control_linkage']
                    ),
                ]),
            ),
            className='fac-Sider'
        )

        # content in right
        self.content = fac.Content(
            className='fac-Content',
            children = fuc.FefferyDiv(
                className='fuc-Grid-container-overview',
                id='DIV_grid_container-overview',
                children = [
                    fuc.FefferyGrid(
                        id = 'FUCGRID_content-overview',
                        className='fuc-Grid',
                        children=[
                            self._init_plot_panels[i].grid_item for i in range( self._n_init_plot_panels )
                        ],
                        layouts= [ 
                            dict(
                                i=p._index, x = (i%3)*16, y = i//3, 
                                w=16, h=1, maxH=1
                            )
                            for i, p in enumerate(self._init_plot_panels)
                        ],
                        cols=48,
                        rowHeight=self._rowHeight_plot_panel,
                        margin=[5, 5],
                        containerPadding=[0,0],
                        autoSize=True,
                        isBounded=True,
                        compactType = 'horizontal',
                    ),
                ]
            ),
        )

        # whole tab
        self.tab = dbc.Tab(
            label = 'Overview',
            tab_id = 'TAB-overview',
            children = [
                self._store,
                fac.Layout(
                    [
                        self.sider,
                        self.content,
                    ],
                    className = 'fac-Layout'
                ),
            ]
        )

clientside_callback(
    ClientsideFunction(
        namespace='overview',
        function_name='grid_container_resize',
    ),
    Output('DIV_grid_container-overview', '_height'),
    Input('DIV_grid_container-overview', '_width'),
    Input('DIV_grid_container-overview', '_height')
)