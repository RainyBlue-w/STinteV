import feffery_antd_components.alias as fac
import dash_mantine_components as dmc

class SelectWithColor:
    
    def __init__(self, index):
        self._index = index
    
        self.select = dmc.Grid(
            gutter=5,
            children=[
                dmc.GridCol(
                    span='auto',
                    children=[
                        fac.Select(
                            locale='en-us',
                            allowClear=False,
                            optionFilterMode='case-sensitive',
                            id={'type': 'SelectWithColor_select', 'index': self._index},
                            style={'width': '100%'}
                        )
                    ]
                ),
                # dmc.GridCol(
                #     span='content',
                #     children=[
                #         fac.ColorPicker(
                #             id={'type': 'SelectWithColor_color', 'index': self._index}
                #         )
                #     ]
                # )
            ]
        )