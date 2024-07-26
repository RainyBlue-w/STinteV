from dash_extensions.enrich import callback, Input, Output, State, no_update, clientside_callback
from dash.exceptions import PreventUpdate
import dash_mantine_components as dmc
import feffery_antd_components.alias as fac
import feffery_utils_components as fuc
from dash_iconify import DashIconify
from dash import html, dcc
from flask_login import current_user

import os

from stintev.config import PathConfig

class UploadPanel:
    
    @staticmethod
    def information_card():
        return dmc.Stack(
            [
                dmc.Tabs(
                    id = 'TABS_upload-dataset',
                    allowTabDeactivation=False,
                    value='new',
                    children=[
                        dmc.TabsList([
                            dmc.TabsTab('New dataset', value='new'),
                            dmc.TabsTab('Existing dataset', value='existing')
                        ]),
                        dmc.TabsPanel(
                            dmc.Stack([
                                dmc.TextInput(
                                    id = 'INPUT_name_new_upload-dataset',
                                    label='Dataset name:',
                                    placeholder='Name of the new dataset',
                                    leftSection=DashIconify(icon='iconoir:database-restore', width=24),
                                ),
                                dmc.TextInput(
                                    id = 'INPUT_description_upload-dataset',
                                    label='Description of :',
                                    placeholder='Description of the dataset',
                                    leftSection=DashIconify(icon='iconoir:page-edit', width=24),
                                ),
                            ]),
                            value='new'
                        ),
                        dmc.TabsPanel(
                            children=dmc.Stack(
                                [
                                    dmc.Select(
                                        id = 'SELECT_name_existing_upload-dataset',
                                        label = 'Select dataset:',
                                        placeholder='Select one from existing private datasets',
                                        leftSection=DashIconify(icon='iconoir:database-star', width=24),
                                        data = sorted(
                                            os.listdir(
                                                os.path.join(PathConfig.DATA_PATH, 'datasets','private',current_user.username)
                                            )
                                        )
                                    )
                                ],
                                gap=0
                            ),
                            value='existing',
                        )
                    ]
                ),
                dmc.Button(
                    children='Continue',
                    rightSection=DashIconify(icon='fluent:next-24-regular', width=24),
                    id='BUTTON_continute_upload-dataset',
                )
            ]
        )

    @staticmethod
    def alert():
        return dmc.Alert(
            'You are browsing as a guest, sign in first to use this!',
            title = 'Permission Denied',
            color = 'violet',
            withCloseButton=False,
            icon = DashIconify(icon='mdi:denied', width=24),
            className='dmc-Alert-guest-dataset'
        )

    @staticmethod
    def upload_card(dataset_name):
        return dmc.Stack(
            [
                fac.DraggerUpload(
                    uploadId=dataset_name,
                    id='UPLOAD_dataset-dataset',
                    locale='en-us',
                    apiUrl='/upload/',
                    fileMaxSize=3072,
                    multiple=False,
                    directory=False,
                    fileTypes = ['h5ad', 'rds'],
                    text='Upload .h5ad or .rds files',
                    hint='Click or drag the files to this area to upload',
                ),
                dmc.Button(
                    children='Finished',
                    color='green',
                    rightSection = DashIconify(icon='fluent:checkmark-circle-24-regular', width=24),
                    id='BUTTON_finished_upload-dataset',
                )
            ]
        )

    def __init__(self):

        self.timeline = dmc.Timeline(
            id='TIMELINE_upload-dataset',
            active=0,
            bulletSize=40,
            lineWidth=2,
            color='green',
            children=[
                dmc.TimelineItem(
                    bullet = DashIconify(icon='solar:document-text-linear', width=32),
                    style={'font-size': '20px'},
                    title = 'Information',
                    children=dmc.Stack(
                        gap=0,
                        className='div-timeline-note-upload',
                        children=[
                            dmc.Text('Create a dataset with basic'),
                            dmc.Text('information')
                        ], 
                    )
                ),
                dmc.TimelineItem(
                    bullet = DashIconify(icon='solar:cloud-upload-linear', width=32),
                    style={'font-size': '20px'},
                    title = 'Upload',
                    children=dmc.Stack(
                        gap=0,
                        className='div-timeline-note-upload',
                        children=[
                            dmc.Text('Upload files'),
                        ], 
                    )
                ),
                # dmc.TimelineItem(
                #     style={'font-size': '20px'},
                #     title = 'Settings',
                #     children=dmc.Stack(
                #         gap=0,
                #         className='div-timeline-note-upload',
                #         children=[
                #             dmc.Text('Configure some settings'),
                #             dmc.Text('for visualization')
                #         ], 
                #     )
                # )
            ]
        )

        self.convert_dialog = html.Div(
            [
                dmc.Modal(
                    closeOnClickOutside = False,
                    size = "60%",
                    title="RDS to H5AD Conversion",
                    id="DIALOG_convert_rds",
                    centered=True,
                    zIndex=10000,
                    children=[
                        dmc.Stepper(
                            id="STEPPER_convert_rds",
                            active=0,
                            children=[
                                dmc.StepperStep(
                                    id="STEP_parse_rds",
                                    loading=False,
                                    label="Parse RDS",
                                    description="Identify file content",
                                    children=[
                                        dmc.Stack(
                                            gap=10,
                                            align="center",
                                            children=[
                                                dmc.Text("", c="red", id='STEP_parse_rds_result'),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl")
                                            ],
                                        )
                                    ]
                                ),
                                dmc.StepperStep(
                                    label="Select Metadata",
                                    description="Filter metadata of interest",
                                    children=[
                                        dmc.Stack(
                                            gap=10,
                                            align="center",
                                            children=[
                                                # dmc.Text("Selecting the metadata you are interested in.", c="#006400"),
                                                dmc.CheckboxGroup(
                                                    id="CHECKBOX_checked_metadata",
                                                    withAsterisk=True,
                                                    mb=10,
                                                    children=[],
                                                    value=[],
                                                ),
                                                dmc.Button(
                                                    "Submit",
                                                    id="BUTTON_submit_metadata",        
                                                    rightSection=DashIconify(icon="formkit:submit"),
                                                ),
                                            ]
                                        )
                                    ]
                                ),
                                dmc.StepperStep(
                                    id="STEP_convert_rds",
                                    loading=False,
                                    label="Convert to H5AD",
                                    description="Generate h5ad file",
                                    children=[
                                        dmc.Stack(
                                            gap=10,
                                            align="center",
                                            children=[
                                                # dmc.Text("Step 3: Converting, please wait......", c="#006400"),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl"),
                                                dmc.Skeleton(h=8, radius="xl")
                                            ],
                                        )
                                    ]
                                ),
                            ],
                        ),
                    ],
                )
            ]
        )

        self.panel = dmc.Group(
            className = 'dmc-Group-upload-dataset',
            children=[
                dcc.Store(id='STORE_dataset_name_upload-dataset'),
                dcc.Store(id='STORE_dataset_description_upload-dataset'),
                self.convert_dialog,
                self.timeline,
                dmc.Card(
                    withBorder=True,
                    shadow="md",
                    radius='lg',
                    className='dmc-card-dataset-upload',
                    id = 'CARD_upload-dataset',
                    children=[UploadPanel.information_card()],
                ) 
            ] 
        ) if current_user.is_authenticated else UploadPanel.alert()

@callback(
    Output('CARD_upload-dataset', 'children'),
    Output('TIMELINE_upload-dataset', 'active'),
    Output('INPUT_name_new_upload-dataset', 'error'),
    Output('INPUT_description_upload-dataset', 'error'),
    Output('SELECT_name_existing_upload-dataset', 'error'),
    Output('STORE_dataset_name_upload-dataset', 'data'),
    Output('STORE_dataset_description_upload-dataset', 'data'),
    
    Input('BUTTON_continute_upload-dataset', 'n_clicks'),
    State('INPUT_name_new_upload-dataset', 'value'),
    State('INPUT_description_upload-dataset', 'value'),
    State('SELECT_name_existing_upload-dataset', 'value'),
    State('TABS_upload-dataset', 'value'),
    prevent_initial_call=True
)
def continute_to_upload(continute, dataset_new, description, dataset_exsiting, active_tab):
    
    if continute:
        if active_tab == 'new':
            if all([dataset_new, description]):
                # 创建dataset目录
                dir_dataset = os.path.join(
                        PathConfig.DATA_PATH, 'datasets', 'private', current_user.username, dataset_new
                    )
                try:
                    os.makedirs(dir_dataset)
                except OSError:
                    pass
                
                # description写入dataset目录下description.txt文件
                f_des = open(os.path.join(dir_dataset, 'description.txt'), 'w')
                f_des.write(description)
                
                return (
                    UploadPanel.upload_card(dataset_new),
                    1, no_update, no_update, no_update,
                    dataset_new, description,
                )
            else:
                return (
                    no_update, no_update, 
                    False if dataset_new else 'Enter dataset name',
                    False if description else 'Enter description',
                    no_update, 
                    no_update, no_update
                )
        elif active_tab == 'existing':
            if dataset_exsiting:
                return (
                    UploadPanel.upload_card(dataset_exsiting),
                    1, no_update, no_update, no_update,
                    dataset_exsiting, no_update,
                )
            else:
                return (
                    no_update, no_update,
                    no_update, no_update,
                    False if dataset_exsiting else 'Select a dataset',
                    no_update, no_update
                )
        else:
            raise PreventUpdate
            
    else:
        raise PreventUpdate

@callback(
    Output('CARD_upload-dataset', 'children'),
    Output('TIMELINE_upload-dataset', 'active'),
    Output('STORE_dataset_name_upload-dataset', 'data'),
    Output('STORE_dataset_description_upload-dataset', 'data'),
    
    Input('BUTTON_finished_upload-dataset', 'n_clicks'),
    prevent_initial_call=True
)
def finished_upload(finished):
    
    if finished:
        return (
            UploadPanel.information_card(), 
            0, None, None
        )
    else:
        raise PreventUpdate