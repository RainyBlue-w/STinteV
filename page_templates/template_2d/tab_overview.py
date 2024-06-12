from dash.dash_table.Format import Format, Group, Scheme, Symbol
from dash_extensions.enrich import html
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac
import feffery_utils_components as fuc

from dash_extensions.enrich import Output, Input, State, callback, clientside_callback, ClientsideFunction, callback_context
from dash import dcc
from dash import ALL, MATCH, Patch, ctx
from dash.exceptions import PreventUpdate

import uuid

from page_templates.template_2d.components.plot_panel import PlotPanel

class TabOverview:
    
    # const
    _width_sider = 250
    _width_sider_collapsed = 60
    _width_drawer = 250
    _rowHeight_plot_panel = 350
    
    # widgets
    
    control_plot_settings = None
    
    control_plot_panels = None
    drawer_plot_panels = None
    _n_init_plot_panels = 2
    _init_plot_panels = [ PlotPanel(uuid.uuid1().hex) for i in range(_n_init_plot_panels) ]
    
    sider = None
    content = None
    tab = None
    
    # Tab init
    def __init__(self) -> None:
        
        self.drawer_plot_panels = html.Div([
            dmc.Drawer(
                id = 'DRAWER_setting_panels-overview-2d',
                className = 'dmc-Drawer-plotPanels',
                size = f'{self._width_drawer}px',
                closeOnClickOutside = False,
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
                        id = 'DRAWER_setting_panels_div-overview-2d',
                        children=[
                            self._init_plot_panels[i].settings for i in range( self._n_init_plot_panels ) 
                        ]
                    ),
                    dmc.Button(
                        "Add panel", id='BUTTON_setting_panels_add-overview-2d',
                        color='teal', fullWidth=True,
                        leftSection=DashIconify(icon="fluent:add-square-20-regular", width=20)
                    )
                ]
            )
        ])
        
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
                        ])
                    ]
                )
            ],
            value='control_plot_settings'
        )
        
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
                                            id = 'NUMBERINPUT_setting_panels-overview-2d',
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
                                        dmc.Button('Setting panels', fullWidth=True, id='BUTTON_setting_panels-overview-2d'),
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

        self.sider = fac.Sider(
            collapsible = False,
            collapsedWidth = self._width_sider_collapsed,
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

        self.content = fac.Content(
            [
                dcc.Store(
                    id='STORE_plotPanelsCurUUID-overview-2d', 
                    data = [
                        self._init_plot_panels[i]._index for i in range( self._n_init_plot_panels )
                    ]
                ),
                fuc.FefferyGrid(
                    id = 'FUCGRID_content-overview-2d',
                    children=[
                        self._init_plot_panels[i].grid_item for i in range( self._n_init_plot_panels )
                    ],
                    layouts= [ 
                        dict(
                            i=p._index, x = (i%3)*2, y = i//3, 
                            w=2, h=1, maxH=1
                        ) 
                        for i, p in enumerate(self._init_plot_panels)
                    ],
                    cols=6,
                    rowHeight=self._rowHeight_plot_panel,
                    margin=[5, 5],
                    containerPadding=[0,0],
                    autoSize=True,
                    isBounded=True,
                    className='fuc-Grid',
                    compactType = 'horizontal',
                )
            ],
            className='fac-Content'
        )

        self.tab = dbc.Tab(
            label = 'Overview',
            tab_id = 'TAB-overview-2d',
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

#region ----Callbacks----

# open the drawer for setting panels
clientside_callback(
    """
    function(n_clicks){
        return true
    }
    """,
    Output('DRAWER_setting_panels-overview-2d', 'opened'),
    Input('BUTTON_setting_panels-overview-2d', 'n_clicks'),
)

# update the rowHeight for PlotPanel-items
@callback(
    output = dict(
        styles = Output({'type': 'PLOTPANEL_item_graph-2d', 'index': ALL}, 'style'),
        rowHeight = Output('FUCGRID_content-overview-2d', 'rowHeight')
    ),
    inputs = [
        Input('NUMBERINPUT_setting_panels-overview-2d', 'value')
    ],
    prevent_initial_call = False,
    suppress_callback_exceptions=True
)
def update_height_for_plot_panel_items(height):

    n_outputs = len(callback_context.outputs_grouping['styles'])
    return {
        'styles': [ {'height': [f'{height-50}px']} ] * n_outputs,
        'rowHeight': height
    }

# add & delete PlotPanel

@callback(
   Output('FUCGRID_content-overview-2d', 'children', allow_duplicate=True),
   Output('DRAWER_setting_panels_div-overview-2d', 'children', allow_duplicate=True),
   Output('STORE_plotPanelsCurUUID-overview-2d', 'data', allow_duplicate=True),
   Output('FUCGRID_content-overview-2d', 'layouts'),
   Output('BUTTON_setting_panels_add-overview-2d', 'disabled'),
   
   Input('BUTTON_setting_panels_add-overview-2d', 'n_clicks'),
   Input({'type': 'PLOTPANEL_settings_button_delete-2d', 'index': ALL}, 'n_clicks'),
   State('STORE_plotPanelsCurUUID-overview-2d', 'data'),
   suppress_callback_exceptions=True,
   prevent_initial_call=True
)
def add_plot_panel(add, delete, uuid_list):
    
    import uuid
    
    tid = ctx.triggered_id
    children_grid = Patch()
    children_drawer = Patch()
    
    if tid == 'BUTTON_setting_panels_add-overview-2d':
        if add and (len(uuid_list) <= 5):
            next_PlotPanel = PlotPanel(uuid.uuid1().hex)
            children_grid.append(
                next_PlotPanel.grid_item
            )
            children_drawer.append(
                next_PlotPanel.settings
            )
            uuid_list.append(
                next_PlotPanel._index
            )
            
    elif ('type' in tid) and (tid['type'] == 'PLOTPANEL_settings_button_delete-2d') and (len(uuid_list) > 1):

        del_index = tid['index']
        for i, index in enumerate(uuid_list):
            if index == del_index:
                del children_grid[i]
                del children_drawer[i]
                del uuid_list[i]
    else:
        raise PreventUpdate

    layouts = []
    for i, index in enumerate(uuid_list):
        layouts.append(
            dict(
                i=index, x = (i%3)*2, y = i//3, 
                w=2, h=1, maxH=1
            ) 
        )

    if len(uuid_list) == 6:
        ban_button = True
    else:
        ban_button = False
        
    return  children_grid, children_drawer, uuid_list, layouts, ban_button