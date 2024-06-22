import dash
import dash._dash_renderer
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_extensions.enrich import html, DashProxy, LogTransform, ServersideOutputTransform, MultiplexerTransform
from dash_extensions.enrich import RedisBackend

from fastapi import FastAPI
from fastapi.middleware.wsgi import WSGIMiddleware
   
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




if __name__ == "__main__":
    
    dashapp.run(
        host='0.0.0.0',
        port=8055,
        debug=True
    )
    
    # import uvicorn
    # app = FastAPI(root_path='/stintev')
    # app.mount("/", WSGIMiddleware(dashapp.server))
    # uvicorn.run(app, port=8055, host='0.0.0.0', log_level='info', debug=True)
