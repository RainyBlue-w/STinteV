import dash_bootstrap_components as dbc
from dash_extensions.enrich import html

from .tab_overview import TabOverview
from .tab_dataset import TabDataset
from .tab_ligand_receptor import TabLigandReceptor
from stintev.config import PathConfig

def render_content(
    path_data_folder: str = PathConfig.DATA_PATH,
):
    
    _tab_overview = TabOverview(path_data_folder)
    _tab_dataset = TabDataset(path_data_folder)
    _tab_ligand_receptor = TabLigandReceptor()

    _tabs = dbc.Tabs(
        children=[
            _tab_dataset.tab,
            _tab_overview.tab,
            _tab_ligand_receptor.tab,
        ],
        active_tab = 'TAB-dataset',
        # active_tab = 'TAB-overview',
        id = 'TABS-all',
    )
    layout = html.Div(
        [
            html.Div(id="notifications-container-main"),
            _tabs,
        ]
    )
    
    return layout