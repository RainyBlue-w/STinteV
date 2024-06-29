from dash import register_page
from dash_extensions.enrich import html

register_page(__name__, path='/login')

layout = html.Div(
    'login'
)