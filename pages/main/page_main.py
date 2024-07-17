import dash_bootstrap_components as dbc
from dash_extensions.enrich import html

from .tab_overview import TabOverview
from .tab_dataset import TabDataset
from stintev.config import PathConfig

import stintev.callbacks.main_c # callbacks

def render_content(
    path_data_folder: str = PathConfig.DATA_PATH,
):
    
    _tab_overview = TabOverview(path_data_folder)
    _tab_dataset = TabDataset(path_data_folder)

    _tabs = dbc.Tabs(
        children=[
            _tab_dataset.tab,
            _tab_overview.tab
        ],
        # active_tab = 'TAB-dataset',
        active_tab = 'TAB-overview',
        id = '_tabs-2d',
    )
    layout = html.Div(
        _tabs
    )
    
    return layout