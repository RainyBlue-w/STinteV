from dash_extensions.enrich import Input, Output, State, no_update
from dash.exceptions import PreventUpdate
from dash import ctx
from flask_login import login_user

import shortuuid

from flask_login import current_user, logout_user, login_user
from stintev import pages
from stintev.components import Modals, Notifications, Navbar
from stintev.server import dashapp
from stintev.models.auth import User

#region 登录系统
@dashapp.callback( # New session ID (sign up)
    Output('DIV_navbar_item_user-app', 'children'), # 刷新navbar sessionID 状态
    Output('app-mount', 'children'), # 刷新页面内容
    Input('BUTTON_new_session_id-app', 'n_clicks'),
)
def new_session_id(n_clicks):
    
    if n_clicks:
        uid = shortuuid.uuid()
        User.add_one_user(uid)
        current_user = User()
        current_user.id = uid
        login_user(current_user)
        return (
            Navbar.navbar_item_user(uid), 
            pages.main.render_content(),
        )

    raise PreventUpdate

@dashapp.callback( # Quit session (logout)
    Output('DIV_navbar_item_user-app', 'children'), # 刷新navbar sessionID 状态
    Output('app-mount', 'children'), # 刷新页面内容
    Input('BUTTON_logout-app', 'n_clicks'),
    prevent_initial_call=True,
)
def quit_session(n_clicks):
    if n_clicks:
        logout_user()
        return (
            Navbar.navbar_item_user(None), 
            pages.main.render_content(),
        )
    raise PreventUpdate

@dashapp.callback( # Have a Session ID (login)
    Output('modal-container', 'children'), # 弹窗收集session ID 输入
    Input('BUTTON_have_session_id-app', 'n_clicks'),
    prevent_initial_call=True,
)
def have_session_id(n_clicks):
    if n_clicks:
        return Modals.modal_have_session_id()
    raise PreventUpdate

@dashapp.callback( # Have a Session ID (login) - Submit
    Output('DIV_navbar_item_user-app', 'children'), # 刷新navbar sessionID 状态
    Output('app-mount', 'children'), # 刷新页面内容
    Output('notifications-container-main', 'children'), # notification, 通知登录状态
    
    Input({'type': 'MODAL_have_session_id', 'index': 'Modal'}, 'okCounts'),
    State({'type': 'MODAL_have_session_id', 'index': 'Input'}, 'value'), # 获取输入的session ID
    prevent_initial_call=True,
)
def have_session_id_submit(okCounts, session_id):

    if okCounts:
        user = User.query_by_id(session_id)
        if user: # 如果 session ID 存在
            login_user(user[0])
            return (
                Navbar.navbar_item_user(session_id),
                pages.main.render_content(),
                []
            )
        else:
            return (
                Navbar.navbar_item_user(None),
                no_update,
                Notifications.notif_no_session_id()
            )
    raise PreventUpdate

#endregion