import dash_mantine_components as dmc
from dash_extensions.enrich import html
from dash_extensions.enrich import Output, Input, MATCH, callback

import os

class DatasetList:
        
    def __init__(
        self, 
        path_data_folder: str, 
        id_prefix: str, # Public | Private
    ):
        
        self.datasets = {
            dataset : [
                file for file in
                sorted(os.listdir(os.path.join(path_data_folder,dataset)))
                if file.endswith('.h5ad')
            ]
            for dataset in 
            sorted(os.listdir(os.path.join(path_data_folder)))
        }
        
        self.list = dmc.CheckboxGroup(
            id={'type': 'DatasetList-checkboxGroup', 'index': id_prefix},
            value=[],
            children=dmc.Grid(
                children=[
                    dmc.GridCol(span=3,
                        children=
                        dmc.Card(
                            radius='md',
                            shadow='sm',
                            withBorder=True,
                            children=[
                                dmc.Checkbox(
                                    value=dataset,
                                    size='lg',
                                    label=[
                                        dmc.Stack(
                                            gap=0,
                                            children=[
                                                dmc.Group(
                                                    children=[
                                                        dmc.Text(
                                                            dataset, size='lg'
                                                        ),
                                                        dmc.Burger(
                                                            id={'type': 'DatasetList-burger-button', 'index': f'{id_prefix}-{dataset}'},
                                                            opened=False,
                                                            style={'position': 'absolute', 'right': '10px'}
                                                        )
                                                    ],
                                                ),
                                                dmc.Space(h=5),
                                                html.Div(
                                                    hidden=True,
                                                    children=dmc.Stack(
                                                        [dmc.Text(file, size='md', c='gray') for file in files],
                                                        gap=0,
                                                    ),
                                                    id={'type': 'DatasetList-dataset-files', 'index': f'{id_prefix}-{dataset}'}
                                                ),
                                            ],
                                        )
                                    ]
                                ),
                                dmc.Space(h=5),
                                dmc.CardSection(
                                    [
                                        dmc.Text('This is an annotation of this dataset', c='gray', style={'margin-left': '10px'})
                                    ]
                                ),
                                dmc.Space(h=2),
                            ],
                        )
                    ) 
                    for dataset, files in self.datasets.items()
                ]
            )
        )
        
@callback(
    Output({'type': 'DatasetList-dataset-files', 'index': MATCH}, 'hidden'),
    Input({'type': 'DatasetList-burger-button', 'index': MATCH}, 'opened')
)
def toggle_hidden(opened):
    return not opened