import dash
import dash._dash_renderer
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_extensions.enrich import html, DashProxy, LogTransform, ServersideOutputTransform, MultiplexerTransform


dash._dash_renderer._set_react_version('18.2.0') # needed for dash_mantine_components v0.14

stylesheets = [
    dbc.themes.BOOTSTRAP,
    "https://unpkg.com/@mantine/dates@7/styles.css",
    "https://unpkg.com/@mantine/code-highlight@7/styles.css",
    "https://unpkg.com/@mantine/charts@7/styles.css",
    "https://unpkg.com/@mantine/carousel@7/styles.css",
    "https://unpkg.com/@mantine/notifications@7/styles.css",
    "https://unpkg.com/@mantine/nprogress@7/styles.css",
]

dash_app = DashProxy(
  __name__, 
  external_stylesheets=stylesheets,
  external_scripts = [
    {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
  ],
  transforms=[
    LogTransform(), ServersideOutputTransform(), MultiplexerTransform()
  ],
  use_pages=True,
)

header = dbc.NavbarSimple(
    [
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem('Test', href='/'),
            ],
            nav=True,
            in_navbar=True,
            label="Dataset",
        ),
    ],
    brand="STinteV",
    color="dark",
    dark=True,
    # sticky='top',
    className='dbc-Navbar-main'
)

dash_app.layout = dmc.MantineProvider(
    [
        dmc.NotificationProvider(),
        header,
        dash.page_container
    ]
)

if __name__ == "__main__":
    dash_app.run(

        host='10.86.60.31',
        port='8055',
        threaded=True,
        proxy=None,

        debug=True
    )