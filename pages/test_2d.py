from dash import register_page
from page_templates.template_2d import layouts
from page_templates.template_2d.callbacks import *

register_page(__name__, path='/')

path_dataset = 'path_dataset'

layout = layouts.init_layout_2d(
    path_dataset=path_dataset
)