import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify
from flask_login import current_user

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
                            dmc.MenuItem('Have a session ID', id='BUTTON_have_session_id-app'),
                            dmc.MenuItem('New Session ID', id='BUTTON_new_session_id-app'),
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
                        )
                    ),
                    dmc.MenuDropdown(
                        [
                            dmc.MenuItem('Quit session', id='BUTTON_logout-app')
                        ]
                    )
                ]
            )


