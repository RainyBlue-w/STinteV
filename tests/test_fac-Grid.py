import dash
import uuid
from dash import html, Patch
import feffery_antd_components as fac
import feffery_utils_components as fuc
from dash.dependencies import Input, Output, State, ALL

app = dash.Dash(__name__)

COLS = 3

app.layout = html.Div(
    [
        fac.AntdCol([
            fac.AntdSpace(id='json-viewer', addSplitLine=True, wrap=True),
            fac.AntdSpace(
                [
                    fac.AntdButton(
                        '添加子项',
                        type='primary',
                        id='add-new-item'
                    ),
                    fac.AntdButton(
                        '删除选中',
                        danger=True,
                        id='delete-selected-items'
                    )
                ],
                style={
                    'marginBottom': 10
                }
            ),
        ]),
        fuc.FefferyGrid(
            [
                fuc.FefferyGridItem(
                    '初始项，不可删除',
                    key='初始项',
                    style={
                        'height': '100%',
                        'display': 'flex',
                        'justifyContent': 'center',
                        'alignItems': 'center'
                    }
                )
            ],
            id='demo-grid',
            layouts=[
                dict(
                    i='初始项',
                    x=0,
                    y=0,
                    w=1,
                    h=1
                )
            ],
            cols=COLS,
            rowHeight=200,
            placeholderBorderRadius='5px',
            margin=[25, 25],
            style={
                'border': '1px dashed #e1dfdd'
            },
            compactType='horizontal'
        )
    ],
    style={
        'padding': '50px 100px'
    }
)

@app.callback(
    [Output('demo-grid', 'children'),
     Output('demo-grid', 'layouts')],
    Input('add-new-item', 'nClicks'),
    State('demo-grid', 'layouts'),
    prevent_initial_call=True
)
def add_new_item(nClicks, origin_layouts):

    new_item_id = str(uuid.uuid4())

    children_patch = Patch()
    layouts_patch = Patch()

    children_patch.append(
        fuc.FefferyGridItem(
            [
                html.Div(
                    new_item_id,
                    style={
                        'height': '100%',
                        'display': 'flex',
                        'justifyContent': 'center',
                        'alignItems': 'center'
                    }
                ),
                fac.AntdCheckbox(
                    checked=False,
                    id={
                        'type': 'delete-current-item',
                        'index': new_item_id
                    },
                    style={
                        'position': 'absolute',
                        'left': 10,
                        'top': 0
                    }
                )
            ],
            key=new_item_id,
            style={
                'position': 'relative'
            }
        )
    )

    layouts_patch.append(
        dict(
            i=new_item_id,
            x=len(origin_layouts) % COLS,
            y=len(origin_layouts) // COLS,
            w=1,
            h=1
        )
    )

    return [
        children_patch,
        layouts_patch
    ]


@app.callback(
    [Output('demo-grid', 'children', allow_duplicate=True),
     Output('demo-grid', 'layouts', allow_duplicate=True)],
    Input('delete-selected-items', 'nClicks'),
    [
        State(
            {
                'type': 'delete-current-item',
                'index': ALL
            },
            'checked'
        ),
        State('demo-grid', 'children'),
        State('demo-grid', 'layouts')
    ],
    prevent_initial_call=True
)
def delete_selected_items(nClicks, checked_list, origin_children, origin_layouts):

    to_delete_items = [
        item['id']['index']
        for item in dash.ctx.states_list[0]
        if item['value']
    ]

    if to_delete_items:
        return [
            [
                child for child in origin_children
                if child['props']['key'] not in to_delete_items
            ],
            [
                layout for layout in origin_layouts
                if layout['i'] not in to_delete_items
            ]
        ]

    return dash.no_update


@app.callback(
    Output('json-viewer', 'children'),
    Input('demo-grid', 'layouts')
)
def view_json_layout(layouts):
    return [
        fuc.FefferyJsonViewer(data = layout)
        for layout in layouts
    ]
    

if __name__ == '__main__':
    app.run(debug=True, port=8057, host='0.0.0.0')
