from dash import register_page

register_page(__name__, path='/')

path_server_folder = '/rad/share/omics-viewer/stintev'

# set layout
from stintev.page_templates.template_2d.init_layout import init_layout_2d
layout = init_layout_2d(
    path_server_folder=path_server_folder
)
