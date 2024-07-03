import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_utils_components as fuc
from dash import html

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
                    w=300,
                    id = 'INPUT_username_signin-login',
                    label='Username:',
                    placeholder='Your username',
                    leftSection=DashIconify(icon='fluent:person-24-regular', width=24),
                    persistence='local'
                ),
                dmc.PasswordInput(
                    w=300,
                    id = 'INPUT_password_signin-login',
                    label="Password:",
                    placeholder="Your password",
                    leftSection=DashIconify(icon="fluent:lock-closed-24-regular", width=24),
                ),
                dmc.Button(
                    'Sign In',
                    id='BUTTON_signin-login',
                    rightSection = DashIconify(icon='mdi:sign-in', width=24),
                    color='blue',
                    w=300,
                ),
                dmc.Button(
                    'Sign Up',
                    id = 'BUTTON_signup-login',
                    rightSection=DashIconify(icon='fluent:document-signature-24-regular', width=24),
                    w=300, color='violet'
                ),
                # cookies for auto-login (md5 password)
                fuc.FefferyCookie(id='COOKIE_username-login',cookieKey='username-login',defaultValue='username'),
                fuc.FefferyCookie(id='COOKIE_password-login',cookieKey='password-login',defaultValue='passname'),
                html.Div(id='debug-login')
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
                        w=300,
                        id = 'INPUT_username_signup-login',
                        label='Username:',
                        placeholder='Your username',
                        leftSection=DashIconify(icon='fluent:person-24-regular', width=24),
                    ),
                    dmc.TextInput(
                        w=300,
                        id= 'INPUT_email_signup-login',
                        label='Email:',
                        placeholder='Your email',
                        leftSection=DashIconify(icon='fluent:mail-24-regular', width=24),
                    ),
                    dmc.PasswordInput(
                        w=300,
                        id = 'INPUT_password_signup-login',
                        label="Password:",
                        placeholder="Your password",
                        leftSection=DashIconify(icon="fluent:lock-closed-24-regular", width=24),
                    ),
                    dmc.PasswordInput(
                        w=300,
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
