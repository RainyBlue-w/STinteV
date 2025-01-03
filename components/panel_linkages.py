from typing import Literal, List
from dash_extensions.enrich import Output, Input, State, no_update, html, callback, ALL, MATCH, ctx
from dash import dcc, Patch
from dash.exceptions import PreventUpdate
from dash_iconify import DashIconify
import dash_mantine_components as dmc

import uuid


class PanelLinkages:

    linkage_types = ['feature', '3D view', 'sample', 'embedding', 'highlighting']
    LinkageTypes = List[Literal['feature', '3D view', 'sample', 'embedding', 'highlighting', None]]
    
    @staticmethod
    def new_linkage(index):
        linkage = dmc.Card(
            style={'width': 'calc(15vw - 30px)'},
            withBorder=True,
            children=[
                dmc.CardSection(
                    children=dmc.Group(
                        [
                            dmc.ThemeIcon(
                                children=DashIconify(icon="mdi:link-variant", width=24),
                                color = 'gray', variant='transparent', size='sm',
                                id = {'type': 'PanelLinkages_icon_linkage', 'index': index},
                            ),
                            dmc.Switch(
                                size='sm', radius='lg', checked=False,
                                onLabel='ON', offLabel='OFF',
                                id={'type': 'PanelLinkages_switch_apply', 'index': index},
                            ),
                            dmc.ActionIcon(
                                DashIconify(icon='gg:close-r'),
                                id = {'type': 'PanelLinkages_button_delete', 'index': index},
                                size='sm', variant='subtle', color='red',
                                style = {'position': 'absolute', 'right': '5px'}
                            ),
                        ],
                        gap = 2,
                        justify='start'
                    ),
                    withBorder=True
                ),
                dmc.MultiSelect(
                    label = 'Items to sync',
                    id = {'type': 'PanelLinkages_select_type', 'index': index},
                    data = [
                        {'label': i, 'value': i}
                        for i in PanelLinkages.linkage_types
                    ],
                ),
                dmc.MultiSelect(
                    label = 'Panels',
                    id = {'type': 'PanelLinkages_select_linkage', 'index': index},
                )
            ]
        )
        return linkage

    @staticmethod
    def linkage_mark(color='gray'):
        return dmc.ThemeIcon(
            children=DashIconify(icon="mdi:link-variant", width=16),
            size = 'sm', color=color, variant='transparent'
        )
    
    @staticmethod
    def linkage_settings():
        init_index = uuid.uuid1().hex
        return dmc.Stack(
            [
                html.Div(
                    children=[PanelLinkages.new_linkage(init_index)],
                    id='DIV_panelLinkages_linkagesContainer',
                ),
                dcc.Store(id='STORE_panelLinkages_linkagesCurUUID', data=[init_index]),
                dmc.Stack(
                    [
                        dmc.Button(
                            'New linkage',
                            leftSection=DashIconify(icon='mdi:link-variant-plus', width=24),
                            id='BUTTON_panelLinkages_add_linkage',
                            fullWidth=True
                        )
                    ]
                )
            ]
        )

@callback( # add & delete linkage
    Output('DIV_panelLinkages_linkagesContainer', 'children'),
    Output('STORE_panelLinkages_linkagesCurUUID', 'data'),
    Output('BUTTON_panelLinkages_add_linkage', 'disabled'),
    
    Input('BUTTON_panelLinkages_add_linkage', 'n_clicks'),
    Input({'type': 'PanelLinkages_button_delete', 'index': ALL}, 'n_clicks'),
    State('STORE_panelLinkages_linkagesCurUUID', 'data'),
    prevent_initial_call=True
)
def add_delete_linkage(click_add, click_delete, list_cur_uuid):
    
    tid = ctx.triggered_id
    children = Patch()
    
    if tid == 'BUTTON_panelLinkages_add_linkage': 
        if click_add and len(list_cur_uuid)<=2: # 最多三个
            uid = uuid.uuid1().hex
            children.append( PanelLinkages.new_linkage(index=uid))
            list_cur_uuid.append(uid)
    elif ('type' in tid) and (tid['type'] == 'PanelLinkages_button_delete'):
        del_index = tid['index']
        for i, index in enumerate(list_cur_uuid):
            if index == del_index:
                del children[i]
                del list_cur_uuid[i]
    else:
        raise PreventUpdate
            
    if len(list_cur_uuid) >= 3:
        ban_add_button = True
    else:
        ban_add_button = False
        
    return children, list_cur_uuid, ban_add_button

