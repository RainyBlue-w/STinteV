import dash_bootstrap_components as dbc
from dash_extensions.enrich import html
from page_templates.template_2d.tab_overview import TabOverview

_tab_overview = TabOverview()

_tabs = dbc.Tabs(
    [
        _tab_overview.tab, 
    ],
    active_tab = 'TAB-overview-2d',
    id = '_tabs-2d'
)

def init_layout_2d(
    path_dataset: str
):
    
    layout_all = html.Div(
        _tabs
    )
    
    return layout_all