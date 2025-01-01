import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac

class Notifications:
    
    @staticmethod
    def notif_no_session_id():
        return fac.Notification(
            message = "Session Failed",
            description = "Session ID doesn't exists!",
            placement = 'top',
            type = 'error'
        )
        
    @staticmethod
    def notif_new_session_id():
        return  fac.Notification(
            message='Save your session ID',
            description='Please save the session ID properly, it cannot be retrieved.',
            placement='top',
            type='info'
        )