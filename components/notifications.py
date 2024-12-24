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
            title = 'Attention',
            message= 'Please save the session ID properly, it cannot be retrieved.',
            color = 'orange',
            action = 'show',
            position='top-center',
            autoClose= 10000,
            icon = DashIconify(icon = 'line-md:alert-twotone', width=24)
        )