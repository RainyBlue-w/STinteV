import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import feffery_utils_components as fuc
from dash import html
from dash_extensions.enrich import clientside_callback, ClientsideFunction, Output, Input, State

card_signin = dmc.Card(
    withBorder=True,
    shadow="md",
    radius='lg',
    className = 'dmc-card-signin-login',
    children=[
        dmc.Stack(
            align='center',
            children = [
                dmc.TextInput(
                    w=350,
                    id = 'INPUT_username_signin-login',
                    label='Username:',
                    placeholder='Your username',
                    leftSection=DashIconify(icon='fluent:person-24-regular', width=24),
                    persistence='local'
                ),
                dmc.PasswordInput(
                    w=350,
                    id = 'INPUT_password_signin-login',
                    label = dmc.Group(
                        gap=0,
                        style={'width': '100%'},
                        children=[
                            'Password:',
                            html.A(
                                children='Forget password?',
                                href='reset-password',
                                style={'position': 'absolute', 'right': '0px'}
                            ),
                        ],
                    ),
                    placeholder="Your password",
                    leftSection=DashIconify(icon="fluent:lock-closed-24-regular", width=24),
                    styles = {'label': {'width': '100%'}},
                    style = {'position': 'relative'}
                ),
                dmc.Button(
                    'Sign In',
                    id='BUTTON_signin-login',
                    rightSection = DashIconify(icon='mdi:sign-in', width=24),
                    color='blue',
                    w=350,
                ),
                dmc.Button(
                    'Sign Up',
                    id = 'BUTTON_signup-login',
                    rightSection=DashIconify(icon='fluent:document-signature-24-regular', width=24),
                    w=350, color='violet'
                ),
            ]
        )
    ]
)

card_signup = dmc.Card(
    withBorder=True,
    shadow="md",
    radius='lg',
    className = 'dmc-card-signup-login',
    children=[
        dmc.Stack(
            align='center',
            children = [
                dmc.Stack([
                    dmc.TextInput(
                        w=350,
                        id = 'INPUT_username_signup-login',
                        label='Username:',
                        placeholder='Your username',
                        leftSection=DashIconify(icon='fluent:person-24-regular', width=24),
                    ),
                    dmc.TextInput(
                        w=350,
                        id= 'INPUT_email_signup-login',
                        label='Email:',
                        placeholder='Your email',
                        leftSection=DashIconify(icon='fluent:mail-24-regular', width=24),
                    ),
                    dmc.Stack(gap='sm', children=[
                        dmc.PasswordInput(
                            w=350,
                            id = 'INPUT_password_signup-login',
                            label="Password:",
                            placeholder="Your password",
                            leftSection=DashIconify(icon="fluent:lock-closed-24-regular", width=24),
                        ),
                        dmc.Progress(
                          id = 'PROGRESS_password_strength_signup-login', value=0, color='red',
                          w=350  
                        ),
                    ]),
                    dmc.PasswordInput(
                        w=350,
                        id = 'INPUT_password_confirm_signup-login',
                        label="Confirm  password:",
                        placeholder="Your password again",
                        leftSection=DashIconify(icon="fluent:lock-closed-24-regular", width=24),
                    ),
                    dmc.SimpleGrid(
                        style={'align-items': 'end'},
                        cols=2,
                        spacing='md',
                        children=[
                            dmc.TextInput(
                                label='Verification Code:',
                                id='INPUT_verification_code-login'
                            ),
                            fuc.FefferyCaptcha(id='CAPTCHA_signup-login') 
                        ]
                    ),
                    dmc.SimpleGrid(
                        cols=2,
                        spacing='sm',
                        children=[
                            dmc.Button(
                                'Return', color='pink', fullWidth=True,
                                id='BUTTON_return_signup-login',
                                leftSection=DashIconify(icon='lets-icons:refund-back', width=24)
                            ),
                            dmc.Button(
                                'Submit', color='green', fullWidth=True,
                                id='BUTTON_enter_signup-login',
                                rightSection=DashIconify(icon='mdi:file-document-box-check-outline', width=24)
                            ),
                        ]
                    ),
                ])
            ]
        )
    ]
)

clientside_callback(
    ClientsideFunction(
        namespace='login',
        function_name='password_strength_zxcvbn',
    ),
    Output('PROGRESS_password_strength_signup-login', 'value'),
    Output('PROGRESS_password_strength_signup-login', 'color'),
    Input('INPUT_password_signup-login', 'value'),
)