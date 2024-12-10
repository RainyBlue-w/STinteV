from os import name
from dash_extensions.enrich import html
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import feffery_antd_components.alias as fac
import feffery_utils_components as fuc
from dash import dcc, no_update, ctx, set_props, MATCH, ALL, Patch
from dash_extensions.enrich import clientside_callback, Output, Input, ClientsideFunction, callback
from dash.exceptions import PreventUpdate

import plotly.express as px
import pandas as pd
import importlib.resources as res
import anndata

import uuid
from flask_login import current_user

from stintev.utils._plot import *
from stintev.utils._methods import cellchatdb_query, bivariate_spatial_association
from stintev.components import SelectWithColor

