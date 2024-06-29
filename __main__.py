import dash
import dash._dash_renderer
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_extensions.enrich import html, DashProxy, LogTransform, ServersideOutputTransform, MultiplexerTransform

import argparse


dash._dash_renderer._set_react_version('18.2.0') # needed for dash_mantine_components v0.14

_stylesheets = [
    dbc.themes.BOOTSTRAP,
    "https://unpkg.com/@mantine/dates@7/styles.css",
    "https://unpkg.com/@mantine/code-highlight@7/styles.css",
    "https://unpkg.com/@mantine/charts@7/styles.css",
    "https://unpkg.com/@mantine/carousel@7/styles.css",
    "https://unpkg.com/@mantine/notifications@7/styles.css",
    "https://unpkg.com/@mantine/nprogress@7/styles.css",
]

_navbar = dbc.NavbarSimple(
    [
        dbc.NavItem(dbc.NavLink('Login', href='login'))
    ],
    brand="STinteV",
    color="dark",
    dark=True,
    className='dbc-Navbar-main'
)

dashapp = DashProxy(
  __name__, 
  external_stylesheets=_stylesheets,
  external_scripts = [
    {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
  ],
  transforms = [
    LogTransform(), ServersideOutputTransform(), MultiplexerTransform()
  ],
  prevent_initial_callbacks=True,
  use_pages=True,
  requests_pathname_prefix='/',
)

dashapp.layout = dmc.MantineProvider(
    [
        dmc.NotificationProvider(),
        _navbar,
        dash.page_container
    ]
)

dashapp.config.suppress_callback_exceptions = True

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run the 3Dviewer app')
    parser.add_argument('--host', default='0.0.0.0', required=False, help='Host to run the app on')
    parser.add_argument('--port', default='8000', required=False, help='Port to run the app on')
    parser.add_argument('--debug', action='store_true', required=False, help='Run app on debug mode')
    args = parser.parse_args()
    
    dashapp.run(
        host=args.host,
        port=args.port,
        debug=args.debug
    )
