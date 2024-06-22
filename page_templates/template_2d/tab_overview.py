from dash.dash_table.Format import Format, Group, Scheme, Symbol
from dash_extensions.enrich import html
import dash_extensions.enrich
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac
import feffery_utils_components as fuc

from dash import dcc

import uuid
import os

from stintev.page_templates._io.read import read_dataset
from stintev.page_templates._plot import *
from stintev.page_templates.template_2d.plot_panel import PlotPanel

class TabOverview():
    
    # const
    _width_sider = 250
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
    def __init__(self, path_server_folder: str) -> None:
                    
        self._init_plot_panels = [ 
            PlotPanel(index=uuid.uuid1().hex),
            PlotPanel(index=uuid.uuid1().hex),
        ]

        # drawer for setting plot panels
        self.drawer_plot_panels = html.Div([
            dmc.Drawer(
                id = 'DRAWER_setting_panels-overview',
                className = 'dmc-Drawer-plotPanels',
                size = f'{self._width_drawer}px',
                closeOnClickOutside = True,
                lockScroll = False,
                withOverlay = False,
                title = dmc.Stack(
                    [
                        dmc.Text('Setting panels', className='dmc-Text-drawerTitle'),
                        dmc.Text('tip: setting the contents of plot panels', className='dmc-Text-drawerSubTitle'),
                    ]
                ),
                children = [ 
                    html.Div(
                        id = 'DRAWER_setting_panels_div-overview',
                        children=[
                            self._init_plot_panels[i].sider_settings for i in range( self._n_init_plot_panels ) 
                        ]
                    ),
                    dmc.Button(
                        "Add panel", id='BUTTON_setting_panels_add-overview',
                        color='teal', fullWidth=True,
                        leftSection=DashIconify(icon="fluent:add-square-20-regular", width=20)
                    )
                ]
            )
        ])

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
                                        dmc.Text('Point size:', className='dmc-Text-label'),
                                        span=5
                                    ),
                                    dmc.GridCol(
                                        dmc.NumberInput(
                                            id='NUMBERINPUT_scatter3dPointsize_3D',
                                            value=3, step=0.5, min=0.1,
                                            persistence=True, persistence_type='local'
                                        ),
                                        span=7
                                    ),
                                ],
                                justify='center', gutter=3, className='dmc-Grid-center'
                            ),
                            dmc.Space(h=5),        
                        ]),
                        # ----Download----
                        html.Div([
                            dmc.Divider(label='Download', labelPosition='center', variant='dashed', className='dmc-divider-sidebar-inline'),
                            dmc.Text('tip: replot to take effect', className='dmc-Text-sidebar-tips'),
                            dmc.Grid(
                                [
                                    dmc.GridCol(
                                        dmc.Select(
                                            label = 'type', id='NUMBERINPUT_scatter3dFigtype_3D',
                                            value='png', data = ['svg', 'png', 'jpeg', 'webp'],
                                            persistence = True, persistence_type = 'local', 
                                        ),
                                        span=6
                                    ),
                                    dmc.GridCol(
                                        dmc.NumberInput(
                                            label = 'resolution', id='NUMBERINPUT_scatter3dFigscale_3D',
                                            value=3, step=1, min=1, 
                                            leftSection=DashIconify(icon='uim:multiply', width=16),
                                            persistence = True, persistence_type = 'local', 
                                        ),
                                        span=6
                                    ),
                                ],
                                justify='center', gutter=3, className='dmc-Grid-center'
                            ),     
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
                                        dmc.Button('Setting panels', fullWidth=True, id='BUTTON_setting_panels-overview'),
                                    ],
                                    span=12
                                )
                            ]
                        ),
                        self.drawer_plot_panels
                    ]
                ),
            ],
            value='control_plot_panels'
        )

        # sider in left
        self.sider = fac.Sider(
            collapsible = False,
            width = self._width_sider,
            children = [
                dmc.Accordion(
                    multiple=True,
                    chevronPosition='right',
                    variant='seperated',
                    children=[
                        self.control_plot_settings,
                        self.control_plot_panels,
                    ]
                )
            ],
            className='fac-Sider'
        )

        # content in right
        self.content = fac.Content(
            [
                dcc.Store(
                    id='STORE_plotPanelsCurUUID-overview', 
                    data = [
                        self._init_plot_panels[i]._index for i in range( self._n_init_plot_panels )
                    ]
                ),
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
                )
            ],
            className='fac-Content'
        )

        # whole tab
        self.tab = dbc.Tab(
            label = 'Overview',
            tab_id = 'TAB-overview',
            children = [
                fac.Layout(
                    [
                        fac.Affix(
                            self.sider
                        ),
                        self.content
                    ],
                    className = 'fac-Layout'
                )
            ]    
        )
