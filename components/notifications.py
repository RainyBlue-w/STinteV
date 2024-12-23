import dash_mantine_components as dmc
from dash_iconify import DashIconify

class Notifications:
    
    @staticmethod
    def notif_no_session_id():
        return  dmc.Notification(
        title='Session Failed',
        message="Session ID doesn't exists!",
        color='red',
        action='show',
        position='top-right',
        autoClose=2000, # ms
        icon = DashIconify(icon='icon-park-outline:close-one', width=24)
    )
        
    @staticmethod
    def notif_new_session_id():
        return  dmc.Notification(
            title = 'New Session ID',
            message= 'The session ID cannot be retrieved!! Please save it properly.',
            color = 'red',
            action = 'show',
            position='top-center',
            autoClose=  False,
            icon = DashIconify(icon = 'iconamoon:attention-circle', width=24)
        )