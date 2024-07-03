import dash._dash_renderer
from dash import dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_extensions.enrich import Input, Output, State, no_update
from dash import html
from dash.exceptions import PreventUpdate
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac

from flask_login import current_user, logout_user
from stintev import pages
from stintev.config import RouterConfig

from stintev.server import dashapp

_navbar = dbc.NavbarSimple(
    [
        dbc.NavItem(id = 'DIV_navbar_item_user-app')
    ],
    brand="STinteV",
    color="dark",
    dark=True,
    className='dbc-Navbar-app'
)

dashapp.layout = dmc.MantineProvider(
    [
        # url监听
        dcc.Location(id='url'),
        
        # 消息提示容器
        dmc.NotificationProvider(),
        html.Div(id="notifications-container"),
        
        # 导航栏
        _navbar,
        
        # 页面内容挂载点
        html.Div(id='app-mount'),
        
        # 路由重定向
        html.Div(id='router-redirect-container')
    ]
)

def _navbar_item_user(user_id: str = None):
    
    if user_id is None:
        return dbc.NavLink('Login', href='login'),
    else:
        return dmc.Menu(
            children=[
                dmc.MenuTarget(
                    dmc.Button(
                        current_user.username,
                        color='gray',
                        variant = 'subtle',
                        leftSection=DashIconify(icon='fluent:person-24-regular', width=24),
                        size='md',
                        className='dmc-Button-user',
                    )
                ),
                dmc.MenuDropdown(
                    [
                        dmc.MenuItem(
                            'Log Out', id='BUTTON_logout-app', n_clicks=0
                        )
                    ]
                )
            ]
        )

@dashapp.callback( # 路由重定向和页面渲染
    Output('DIV_navbar_item_user-app', 'children'),
    Output('app-mount', 'children'),
    Output('router-redirect-container', 'children'),
    
    Input('url', 'pathname'),
)
def router(pathname):
    
    # 过滤非法pathname
    if pathname not in RouterConfig.VALID_PATHNAME:
        return (
            _navbar_item_user(None),
            pages._404.render_content(),
            None
        )
    
    # 检查是否已登录
    if current_user.is_authenticated: # 若已登录
        if pathname == '/login': # 重定向回主页面
            return (
                dash.no_update,
                dash.no_update,
                dcc.Location(pathname='/', id='router-redirect'),
            )
        else: # 登录状态的主页面
            return (
                _navbar_item_user(current_user.username),
                pages.main.render_content(),
                None,
            )
    else:  # 若未登录
        if pathname == '/login': # 渲染登录页面
            return (
                _navbar_item_user(None),
                pages.login.render_content(),
                None
            )
        else: # 未登录状态的主页面
            return (
                _navbar_item_user(None),
                pages.main.render_content(),
                dcc.Location(pathname='/', id='router-redirect'),
            )

@dashapp.callback( # 登出
    Output('DIV_navbar_item_user-app', 'children'),
    Output('app-mount', 'children'),
    Output('router-redirect-container', 'children'),

    Input('BUTTON_logout-app', 'n_clicks')
)
def logout_app(n_clicks):
    if n_clicks:
        logout_user()
        return (
            _navbar_item_user(None),
            pages.main.render_content(),
            dcc.Location(pathname='/', id='router-redirect'),
        )
    else:
        raise PreventUpdate