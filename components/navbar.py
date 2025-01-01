import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify
from flask_login import current_user
from dash_extensions.enrich import Output, Input, State, ClientsideFunction, clientside_callback

class Navbar:
    
    @staticmethod
    def navbar():
        _navbar = dbc.NavbarSimple(
            [
                dbc.NavItem(id = 'DIV_navbar_item_user-app')
            ],
            brand="STinteV",
            color="dark",
            dark=True,
            className='dbc-Navbar-app'
        )
        return _navbar

    @staticmethod
    def navbar_item_user(user_id: str | None = None):
        if user_id is None:
            return dmc.Menu(
                children = [
                    dmc.MenuTarget(
                        dmc.Button(
                            'No session ID', 
                            color='gray',
                            variant = 'subtle',
                            leftSection=DashIconify(icon='material-symbols:add-ad', width=24),
                            size='md',
                            className='dmc-Button-user',
                        )
                    ),
                    dmc.MenuDropdown(
                        [
                            dmc.MenuLabel('Session ID'),
                            dmc.MenuItem(
                                'Have a Session ID', id='BUTTON_have_session_id-app',
                                leftSection=DashIconify(icon='hugeicons:id', width=20)
                            ),
                            dmc.MenuItem(
                                'New Session ID', id='BUTTON_new_session_id-app',
                                leftSection=DashIconify(icon='material-symbols:add-ad', width=20)
                            ),
                        ]
                    )
                ]
            )
        else:
            return dmc.Menu(
                children=[
                    dmc.MenuTarget(
                        dmc.Button(
                            f'{current_user.id}',
                            color='gray',
                            variant = 'subtle',
                            leftSection=DashIconify(icon='material-symbols:ad', width=24),
                            size='md',
                            className='dmc-Button-user',
                            id = 'BUTTON_show_session_id-app'
                        )
                    ),
                    dmc.MenuDropdown(
                        dmc.Stack(
                            [
                                dmc.MenuLabel('Session ID'),
                                dmc.Text(f'{current_user.id}', className='navbar-menu-session-id'),
                                dmc.MenuLabel('Operation'),
                                dmc.MenuItem(
                                    'Quit session', id='BUTTON_logout-app',
                                    leftSection=DashIconify(icon='ci:exit', width=20),
                                )
                            ]
                        )
                    )
                ]
            )
