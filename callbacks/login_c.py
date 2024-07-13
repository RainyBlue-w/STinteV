from email import message
from dash_extensions.enrich import Input, Output, State, no_update
from dash.exceptions import PreventUpdate
from dash import dcc
from dash_iconify import DashIconify
from flask_login import login_user
import dash_mantine_components as dmc

from email_validator import validate_email, EmailNotValidError

from stintev.server import dashapp
from stintev.models.auth import UserAccount, User
from stintev.utils import str2md5
from stintev.components import login_card

@dashapp.callback( # Sign In 
    Output('INPUT_username_signin-login', 'error'),
    Output('INPUT_password_signin-login', 'error'),
    Output('login-redirect-container', 'children'),
    
    Input('BUTTON_signin-login', 'n_clicks'),
    
    State('INPUT_username_signin-login', 'value'),
    State('INPUT_password_signin-login', 'value'),
    prevent_initial_call=True,
)
def signin(n_clicks, username, password):
    
    if n_clicks:
        # 输入值全不为空
        if all([n_clicks, username, password]):
            query_user_result = UserAccount.query_user(username)
            if query_user_result: # 如果用户名存在
                if query_user_result[0]['password'] == str2md5(password):
                    # 密码正确
                    current_user = User()
                    current_user.id = username
                    current_user.username = username
                    login_user(current_user)
                    return  (
                        False,
                        False,
                        dcc.Location(pathname='/', id='login-redirect')
                    )
                else:  # 密码错误
                    return False, 'Invalid password', None
            else: # 用户名不存在
                return  'Unknown user name', False, None
        return [
            False if username else 'Enter username',
            False if password else 'Enter password',
            None
        ]
    raise PreventUpdate

@dashapp.callback( # Sign Up
    Output('INPUT_username_signup-login', 'error'),
    Output('INPUT_email_signup-login', 'error'),
    Output('INPUT_password_signup-login', 'error'),
    Output('INPUT_password_confirm_signup-login', 'error'),
    Output('INPUT_verification_code-login', 'error'),
    Output('notifications-container', 'children'), # notification
    Output('CAPTCHA_signup-login', 'refresh'),
    Output('INPUT_verification_code-login', 'value'),
    Output('DIV_card-login', 'children'),
    
    Input('BUTTON_enter_signup-login', 'n_clicks'),
    
    State('INPUT_username_signup-login', 'value'),
    State('INPUT_email_signup-login', 'value'),
    State('INPUT_password_signup-login', 'value'),
    State('PROGRESS_password_strength_signup-login', 'value'),
    State('INPUT_password_confirm_signup-login', 'value'),
    State('INPUT_verification_code-login', 'value'),
    State('CAPTCHA_signup-login', 'captcha')
)
def signup(n_clicks, username, email, password, strength, confirm, veri_code, captcha):
    
    if n_clicks:
        # 输入值全不为空
        if all([n_clicks, username, email, password, confirm, veri_code]):
            
            if UserAccount.query_user(username):
                return 'Username already exists', False, False, False, False, [], True, '', no_update
            
            try:
                validate_email(email, check_deliverability=True)
            except EmailNotValidError as e:
                return False, 'Invalid email', False, False, False, [], True, '', no_update
            
            if strength < 60:
                return False, False, 'Password too weak', False, False, [], True, '', no_update
            if password != confirm:
                return False, False, 'Password mismatch', 'Password mismatch', False, [], True, '', no_update
            if captcha != veri_code:
                return False, False, False, False, 'Invalid verification code', [], True, '', no_update
            UserAccount.add_one_user(username, email, password, 'normal')
            note = dmc.Notification(
                title = 'Registration Successful!',
                action = 'show',
                message = f"You have successfully registered as normal user '{username}'!",
                icon=DashIconify(icon='ic:round-celebration', width=24)
            )
            return False, False, False, False, False, note, True, '', login_card.card_signin
        return [
            False if username else 'Enter username',
            False if email else 'Enter email',
            False if password else 'Enter password',
            False if confirm else 'Enter password again',
            False if veri_code else 'Enter verification code',
            [],
            True,
            '',
            no_update
        ]
    raise PreventUpdate
