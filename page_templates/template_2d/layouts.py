from dash import Dash
from dash import dcc, html, dash_table, no_update, State, Patch
from dash import DiskcacheManager, clientside_callback, ctx, ClientsideFunction
from dash import ALL, MATCH, ALLSMALLER
from dash.dash_table.Format import Format, Group, Scheme, Symbol
from dash.exceptions import PreventUpdate

from dash_extensions.enrich import Output, Input, html, callback, DashProxy
from dash_extensions.enrich import MultiplexerTransform, Trigger, TriggerTransform


import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac
import feffery_utils_components as fuc

import plotly.express as px
import plotly.graph_objects as go

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import math
import os
from functools import reduce

from typing import List, Dict
import diskcache


def init_layout_2d(
    path_dataset: str
):
    
    layout_all = html.Div(
        'Test'
    )
    
    return layout_all