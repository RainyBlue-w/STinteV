import dash._dash_renderer
from dash import dcc
import dash_mantine_components as dmc
from dash_extensions.enrich import Input, Output, State, no_update
from dash import html

from flask_login import current_user
from stintev import pages
from stintev.config import RouterConfig
from stintev.components import Navbar
from stintev.server import dashapp
from stintev.models.auth import User

from stintev.callbacks import login_c, main_c


dashapp.layout = dmc.MantineProvider(
    [
        # url监听
        dcc.Location(id='url'),
        
        # 弹窗modal容器
        html.Div(id='modal-container'),
        
        # 导航栏
        Navbar.navbar(),
        
        # 页面内容挂载点
        html.Div(id='app-mount'),
        
        # 路由重定向
        html.Div(id='router-redirect-container')
    ]
)

@dashapp.callback( # 路由重定向和页面渲染
    Output('DIV_navbar_item_user-app', 'children'),
    Output('app-mount', 'children'),
    Output('router-redirect-container', 'children'),
    
    Input('url', 'pathname'),
    Input('url', 'search')
)
def router(pathname, search):
    
    # print(f'pathname:{pathname}, search:{search}')
    
    # / /login /reset_password
    # 过滤非法pathname
    if pathname not in RouterConfig.VALID_PATHNAME:
        return (
            Navbar.navbar_item_user(None),
            pages._404.render_content(),
            None
        )
    
    if pathname == '/login':
        if current_user.is_authenticated:
            return (
                dash.no_update,
                dash.no_update,
                dcc.Location(pathname='/', id='router-redirect'),
            )
        else:
            return (
                Navbar.navbar_item_user(None),
                pages.login.render_content(),
                dcc.Location(pathname='/login', id='router-redirect'),
            )
    
    if pathname == '/reset-password':
        if search: # reset password
            reset_token = search.split('=')[-1]
            user = User.verify_reset_token(reset_token) 
            if user: # token valid
                return (
                    Navbar.navbar_item_user(None),
                    pages.reset_password.render_content(user_id=user.id),
                    None,
                )
            else: # token invalid
                return (
                    Navbar.navbar_item_user(None),
                    pages._404.render_content(),
                    None,
                )
        else: # request reset
            return (
                Navbar.navbar_item_user(None),
                pages.reset_password.render_content(user_id=None),
                dcc.Location(pathname='/reset-password', id='router-redirect'),
            )
    
    else:
        if current_user.is_authenticated:
            return (
                Navbar.navbar_item_user(current_user.id),
                pages.main.render_content(),                                            
                None,
            )
        else:
            return (
                Navbar.navbar_item_user(None),
                pages.main.render_content(),
                dcc.Location(pathname='/', id='router-redirect'),
            )


