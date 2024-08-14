from dash import html
import feffery_antd_components as fac


def render_content():

    return html.Div(
        [
            html.Div(
                [
                    fac.AntdResult(
                        status='404',
                        title="404 not found",
                        subTitle='Sorry, the page you visited does not exist.',
                        style={
                            'paddingBottom': 0,
                            'paddingTop': 0
                        }
                    ),
                    fac.AntdButton(
                        'Back to main page',
                        type='link',
                        href='/',
                        target='_self'
                    )
                ],
                style={
                    'textAlign': 'center'
                }
            )
        ],
        style={
            'height': '100vh',
            'display': 'flex',
            'alignItems': 'center',
            'justifyContent': 'center'
        }
    )
