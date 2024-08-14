from dash_extensions.enrich import html, Output, Input, State, callback, no_update
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_utils_components as fuc
from stintev.components import ResetpwdCard
from stintev.models.auth import User
import stintev.callbacks.reset_password_c

def render_content(
    user_id: str | None = None,
):
    
    # requeset reset
    if user_id is None:
        layout = html.Div(
            className='div-top-request-reset',
            children = [
                html.Div(id='request-reset-redirect-container'),
                html.Div(
                    ResetpwdCard.request_card(), id='DIV_card-request-reset',
                    style={
                        'width': '100%', 'display': 'flex', 
                        'justify-content': 'center'
                    }
                ),
            ]
        )
    
    # reset password
    else:
        layout = html.Div(
            className='div-top-resetpwd',
            children = [
                html.Div(id='resetpwd-redirect-container'),
                html.Div(
                    ResetpwdCard.reset_card(user_id=user_id), id='DIV_card-resetpwd',
                    style={
                        'width': '100%', 'display': 'flex', 
                        'justify-content': 'center'
                    }
                ),
            ]
        )
    
    return layout