import dash._dash_renderer
import dash_bootstrap_components as dbc
from dash_extensions.enrich import DashProxy, LogTransform, ServersideOutputTransform, MultiplexerTransform

from flask_login import LoginManager, current_user
from flask import request

import os

from stintev.models.auth import User, UserAccount
from stintev.config import PathConfig

dash._dash_renderer._set_react_version('18.2.0') # needed for dash_mantine_components v0.14

_stylesheets = [
    dbc.themes.BOOTSTRAP,
    "https://unpkg.com/@mantine/dates@7/styles.css",
    "https://unpkg.com/@mantine/code-highlight@7/styles.css",
    "https://unpkg.com/@mantine/charts@7/styles.css",
    "https://unpkg.com/@mantine/carousel@7/styles.css",
    "https://unpkg.com/@mantine/notifications@7/styles.css",
    "https://unpkg.com/@mantine/nprogress@7/styles.css",
]

dashapp = DashProxy(
  __name__, 
  external_stylesheets=_stylesheets,
  external_scripts = [
    {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
  ],
  transforms = [
    LogTransform(), ServersideOutputTransform(), MultiplexerTransform()
  ],
  prevent_initial_callbacks=True,
  requests_pathname_prefix='/',
)

dashapp.config.suppress_callback_exceptions = True
dashapp.server.secret_key = '4e17572a360111efaa763cecef387085'

# login

login_manager = LoginManager()
login_manager.login_view = '/'
login_manager.init_app(dashapp.server)

# 创建回调函数，作用是供login_manager通过username信息获取对应的User对象
@login_manager.user_loader
def load_user(username):
    if UserAccount.query_user(username):
        curr_user = User()
        curr_user.id = username
        curr_user.username = username
        return curr_user


# upload

@dashapp.server.route('/upload/', methods=['POST'])
def upload():
    '''
    构建文件上传服务
    :return:
    '''

    # 获取上传id参数，用于指向保存路径
    # uploadId <-> dataset
    uploadId = request.values.get('uploadId')

    # 获取上传的文件名称
    filename = request.files['file'].filename

    # 基于用户名，若本地不存在则会自动创建目录
    path_to_save = os.path.join(PathConfig.DATA_PATH, 'datasets', 'private', current_user.username, uploadId)
    try:
        os.makedirs(path_to_save)
    except OSError:
        pass

    # 流式写出文件到指定目录
    with open(os.path.join(path_to_save, filename), 'wb') as f:
        # 流式写出大型文件，这里的10代表10MB
        for chunk in iter(lambda: request.files['file'].read(1024 * 1024 * 10), b''):
            f.write(chunk)

    return {'status': 'success', 'message': 'upload success', 'filename': filename}