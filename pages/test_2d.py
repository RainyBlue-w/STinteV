from dash import register_page

register_page(__name__, path='/')

path_dataset = 'path_dataset'

# set layout
from page_templates.template_2d.init_layout import init_layout_2d
layout = init_layout_2d(
    path_dataset=path_dataset
)