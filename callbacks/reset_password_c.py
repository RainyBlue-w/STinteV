from typing import Dict

from dash_extensions.enrich import Input, Output, State, no_update
from dash.exceptions import PreventUpdate
from dash import dcc
from dash_iconify import DashIconify
import dash_mantine_components as dmc

from flask_mail import Message

from stintev.server import dashapp, mail
from stintev.models.auth import User
from stintev.config import NetConfig

# request reset
def _send_reset_email(user: User):
  token = user.get_reset_token()
  msg = Message(
    'test subject',
    sender = dashapp.server.config['MAIL_USERNAME'],
    recipients = [user.email],
  )
  msg.body = f'''To reset your password, visit the following link in 10 minutes:
  \n
  {NetConfig.SERVER_ADDRESS}/reset-password?token={token}
  \n
  if you did not make this request then simply ignore this email and no changes will be made.
  '''
  mail.send(msg)
  
@dashapp.callback(
    Output('notifications-container', 'children'),
    
    Input('BUTTON_send-requeset-reset', 'n_clicks'),
    State('INPUT_email-request-reset', 'value'),
    prevent_initial_call=True,
)
def send_reset_email(click_send, email):
    
    if click_send:
        user = User.query_by_email(email)[0]
        note = dmc.Notification(
            title = f"Email sent to {user.email}",
            action = 'show',
            message = f"{user.id}, an email has been sent to {user.email} with instructions to reset your password",
            icon=DashIconify(icon='iconoir:send-mail', width=24)
        )
        _send_reset_email(user)
        return note
    
    raise PreventUpdate
  
# reset password
@dashapp.callback(
	Output('notifications-container', 'children'),
	Output('INPUT_password-resetpwd', 'error'),
	Output('INPUT_password_confirm-resetpwd', 'error'),
	Output('resetpwd-redirect-container', 'children'),
 
	Input('BUTTON_enter-resetpwd', 'n_clicks'),
	State('url', 'search'),
	State('INPUT_password-resetpwd', 'value'),
	State('INPUT_password_confirm-resetpwd', 'value'),
	State('PROGRESS_password_strength-resetpwd', 'value'),
	prevent_initial_call=True,
)
def confirm_reset(click_enter, search, password, password_confirm, strength):
	if click_enter:
     
		if strength < 60:
			return (
				no_update,
				'Password too weak',
				False,
				no_update
			)
   
		if password != password_confirm:
			return (
				no_update,
				False,
				'Password mismatch',
				no_update
			)
   
		reset_token = search.split('=')[-1]
		user = User.verify_reset_token(reset_token)
		update = User.update_password(user_id=user.id, new_password=password)
		if update:
			note = dmc.Notification(
				title = 'Password reset',
				action = 'show',
				message = 'Your password has been reset successfully',
				icon = DashIconify(icon='solar:shield-check-linear', width=24)
			)
			return (
				note,
				False,
				False,
				dcc.Location(pathname='/login', id='router-redirect')
			)
		else:	
			note = dmc.Notification(
				title = 'Password reset failed',
				action = 'show',
				message = 'Your password reset failed',
				icon = DashIconify(icon='solar:shield-cross-linear', width=24)
			)
			return (
				note,
				False,
				False,
				no_update
			)
	
	raise PreventUpdate