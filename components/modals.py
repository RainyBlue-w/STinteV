import feffery_antd_components.alias as fac

class Modals:
    
    @staticmethod
    def modal_have_session_id():
        return  fac.Modal(
            id = {'type': 'MODAL_have_session_id', 'index': 'Modal'},
            title = 'Existing session ID:',
            children = [
                fac.Input(
                    id = {'type': 'MODAL_have_session_id', 'index': 'Input'},
                    placeholder='Enter the session ID you have',
                ),
            ],
            renderFooter=True, okText='Submit', cancelText='Cancel',
            cancelButtonProps={'danger': True},
            visible=True,
        )
