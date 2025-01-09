import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac
from dash import set_props, ctx

class Notifications:
    
    @staticmethod
    def handle_global_error(e):
        set_props(
            'notifications-container',
            {
                'children': []
                # 'children': fac.Message(
                #     content = f'id:{ctx.triggered_id}, error:{e}',
                #     type = 'error'
                # )
            }
        ) 
    
    @staticmethod
    def notif_no_session_id():
        return fac.Notification(
            message = "Session Failed",
            description = "Session ID doesn't exists!",
            placement = 'top',
            duration = 5,
            type = 'error'
        )
        
    @staticmethod
    def notif_new_session_id():
        return  fac.Notification(
            message='Save your session ID',
            description='Please save the session ID properly, it cannot be retrieved.',
            placement='top',
            duration = 5,
            type='info'
        )
        
    @staticmethod
    def notif_metadata_category_exceed():
        return fac.Notification(
            message='Too many categories',
            description='The number of categories in metadata exceeds the maximum limit (100).',
            placement='top',
            duration = 5,
            type='error'
        )