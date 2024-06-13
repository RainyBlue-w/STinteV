import dash_bootstrap_components as dbc
from dash_extensions.enrich import html
from .tab_overview import TabOverview


def init_layout_2d(
    path_dataset: str
):
    
    _tab_overview = TabOverview(path_dataset)

    _tabs = dbc.Tabs(
        [
            _tab_overview.tab, 
        ],
        active_tab = 'TAB-overview-2d',
        id = '_tabs-2d'
    )
    layout_all = html.Div(
        _tabs
    )
    
    return layout_all