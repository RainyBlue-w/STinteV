import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import feffery_utils_components as fuc
from dash import html
from dash_extensions.enrich import clientside_callback, ClientsideFunction, Output, Input, State
from stintev.models.auth import User

class ResetpwdCard:
    
    @staticmethod
    def request_card():
        card_request = dmc.Card(
            withBorder=True,
            shadow="md",
            radius='lg',
            className = 'dmc-card-signup-login',
            children = dmc.Stack(
                children=[
                    dmc.TextInput(
                        label = 'Your Email:',
                        description = 'The email you used to register your account',
                        w = 350,
                        id = 'INPUT_email-request-reset',
                    ),
                    dmc.Button(
                        'Send Email', color='blue', fullWidth=True,
                        id='BUTTON_send-requeset-reset',
                        leftSection=DashIconify(icon='fluent:mail-24-regular', width=24)
                    ),
                ]
            )
        )
        return card_request

    
    @staticmethod
    def reset_card(user_id: str | None = None):
        card_reset = dmc.Card(
            withBorder=True,
            shadow="md",
            radius='lg',
            className = 'dmc-card-signup-login',
            children = [
                dmc.CardSection(
                    dmc.Group(
                        children=[
                            dmc.Text('Reset password', size='md'),
                            dmc.Button(
                                user_id, variant='subtle', color='gray',
                                rightSection=DashIconify(icon='fluent:person-24-regular', width=20),
                            )
                        ],
                        justify='space-between'
                    ),
                    withBorder=True,
                    inheritPadding=True,
                    py='xs',
                ),
                dmc.Space(h=10),
                dmc.Stack(
                    [
                        dmc.Stack(gap='sm', children=[
                            dmc.PasswordInput(
                                w=350,
                                id = 'INPUT_password-resetpwd',
                                label="New password:",
                                placeholder="Your password",
                                leftSection=DashIconify(icon="fluent:lock-closed-24-regular", width=24),
                            ),
                            dmc.Progress(
                                id = 'PROGRESS_password_strength-resetpwd', value=0, color='red',
                                w=350  
                            ),
                        ]),
                        dmc.PasswordInput(
                            w=350,
                            id = 'INPUT_password_confirm-resetpwd',
                            label="Confirm  password:",
                            placeholder="Your password again",
                            leftSection=DashIconify(icon="fluent:lock-closed-24-regular", width=24),
                        ),
                        dmc.Button(
                            'Submit', color='green', fullWidth=True,
                            id='BUTTON_enter-resetpwd',
                            rightSection=DashIconify(icon='mdi:file-document-box-check-outline', width=24)
                        ),
                    ]
                )
            ]
        )
        
        return card_reset
    
clientside_callback(
    ClientsideFunction(
        namespace='login',
        function_name='password_strength_zxcvbn',
    ),
    Output('PROGRESS_password_strength-resetpwd', 'value'),
    Output('PROGRESS_password_strength-resetpwd', 'color'),
    Input('INPUT_password-resetpwd', 'value'),
)