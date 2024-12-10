from dash_extensions.enrich import html, Output, Input, State, callback, no_update
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_utils_components as fuc

from stintev.components import login_card
import stintev.callbacks.login_c # callbacks

def render_content():
    layout = html.Div(
        className='div-top-login',
        children = [
            html.Div(id='login-redirect-container'),
            html.Div(
                login_card.card_signin, id='DIV_card-login',
                style={'width': '100%', 'display': 'flex', 'justify-content': 'center'}
            ),
        ]
    )
    return layout

@callback(
    Output('DIV_card-login', 'children'),
    Input('BUTTON_signup-login', 'n_clicks'),
    prevent_initial_call=True
)
def switch_to_signup(n_clicks):
    if n_clicks:
        return login_card.card_signup
    return no_update

@callback(
    Output('DIV_card-login', 'children'),
    Input('BUTTON_return_signup-login', 'n_clicks'),
    prevent_initial_call=True,
)
def switch_to_signin(n_clicks_return):
    if n_clicks_return:
        return login_card.card_signin
    return no_update
