from dash import register_page

register_page(__name__, path='/')

path_dataset = '/rad/share/omics-viewer/spatial/matrix_data'

# set layout
from stintev.page_templates.template_2d.init_layout import init_layout_2d
layout = init_layout_2d(
    path_dataset=path_dataset
)