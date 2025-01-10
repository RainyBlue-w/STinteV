import dash._dash_renderer
import dash_bootstrap_components as dbc
from dash_extensions.enrich import DashProxy, LogTransform, ServersideOutputTransform
from dash_extensions.enrich import MultiplexerTransform, TriggerTransform, CycleBreakerTransform
from dash_extensions.enrich import RedisCache, FileSystemCache, RedisBackend

from flask import send_from_directory, request
from flask_login import LoginManager, current_user
from flask import request
from flask_mail import Mail
import os
from uuid import uuid4

from stintev.models.auth import User
from stintev.config import *
from stintev.forms import RequesetResetPwdForm, ResetPwdForm
from stintev.components import Notifications

from dash import DiskcacheManager, CeleryManager
import diskcache

launch_uid = uuid4()

if BackendConfig.BACKGROUND_CALLBACK_ON:
    from celery import Celery
    celery_app = Celery(
        __name__, 
        broker=f'redis://:{SecretConfig.REDIS_PASSWORD}@{RedisConfig.HOST}:{RedisConfig.PORT}/0', 
        backend=f'redis://:{SecretConfig.REDIS_PASSWORD}@{RedisConfig.HOST}:{RedisConfig.PORT}/0'
    )
    background_callback_manager = CeleryManager(
        celery_app, 
        # cache_by=[lambda: launch_uid], expire=CeleryConfig.EXPIRE_TIME
    ) # Don't expire for current user for 24 hours
else:
    # Diskcache for non-production apps when developing locally
    import diskcache
    cache = diskcache.Cache(PathConfig.CACHE_PATH)
    background_callback_manager = DiskcacheManager(
        # cache, cache_by=[lambda: launch_uid], expire=60*60
    ) # Don't expire for current user for 1 hours

if BackendConfig.DASH_EXTENSIONS_ON:
    dash_extensions_backend = RedisBackend(
        host=RedisConfig.HOST, 
        port=RedisConfig.PORT, 
        password=SecretConfig.REDIS_PASSWORD
    )
else:
    dash_extensions_backend = FileSystemCache(cache_dir=PathConfig.CACHE_PATH)

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

_external_scripts = [
  'https://cdnjs.cloudflare.com/ajax/libs/cookie-banner/1.2.2/cookiebanner.min.js',
  {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'},
]

dashapp = DashProxy(
  __name__,
  title = 'STinteV',
  external_stylesheets=_stylesheets,
  external_scripts = _external_scripts,
  prevent_initial_callbacks=True,
  requests_pathname_prefix='/',
  suppress_callback_exceptions = True,
  transforms=[MultiplexerTransform(), ServersideOutputTransform(backends=[dash_extensions_backend])],
  background_callback_manager=background_callback_manager,
  on_error = Notifications.handle_global_error
)

dashapp.index_string = '''
<!DOCTYPE html>
<html>
    <script>
    var _hmt = _hmt || [];
    (function() {
    var hm = document.createElement("script");
    hm.src = "https://hm.baidu.com/hm.js?afcaf8eaef1c6a82ef6d1a928719988c";
    var s = document.getElementsByTagName("script")[0]; 
    s.parentNode.insertBefore(hm, s);
    })();
    </script>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
        <script type="text/javascript" id="cookiebanner"
            src="https://cdn.jsdelivr.net/gh/dobarkod/cookie-banner@1.2.2/dist/cookiebanner.min.js"
            data-height="36px" data-position="bottom" data-font-size="18px"
            data-close-text="Got it!"
            data-message="We use cookies to preserve your session state, by using this website you agree to our use of cookies.">
        </script>
    </body>
</html>
'''

# config

dashapp.config.suppress_callback_exceptions = True
dashapp.server.secret_key = SecretConfig.SECRET_KEY
dashapp.server.config['SECRET_KEY'] = SecretConfig.SECRET_KEY

# mail box to send reset password email from
dashapp.server.config['MAIL_SERVER'] = 'smtp.163.com'
dashapp.server.config['MAIL_PORT'] = 465
dashapp.server.config['MAIL_USE_SSL'] = True
dashapp.server.config['MAIL_USERNAME'] = 'stintev@163.com'
dashapp.server.config['MAIL_PASSWORD'] = 'XDMCUELZOKCIKWTY'

# flask-mail
mail = Mail(dashapp.server)

# login

login_manager = LoginManager()
login_manager.login_view = '/'
login_manager.init_app(dashapp.server)

# 创建回调函数，作用是供login_manager通过id信息获取对应的User对象
@login_manager.user_loader
def load_user(user_id):
    if User.query_by_id(user_id):
        curr_user = User()
        curr_user.id = user_id
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
    path_to_save = os.path.join(PathConfig.DATA_PATH, 'datasets', 'private', current_user.id, uploadId)
    try:
        os.makedirs(path_to_save)
    except OSError:
        pass

    # 流式写出文件到指定目录
    with open(os.path.join(path_to_save, filename), 'wb') as f:
        # 流式写出大型文件，这里的100代表100MB
        for chunk in iter(lambda: request.files['file'].read(1024 * 1024 * 50), b''):
            f.write(chunk)

    return {'status': 'success', 'message': 'upload success', 'filename': filename}

# download
@dashapp.server.route('/download_public/', methods=['GET'])
def download_public():
    dataset = request.args.get('dataset')
    file = request.args.get('file')
    return send_from_directory(f'{PathConfig.DATA_PATH}/datasets/public/{dataset}', file)

# download
@dashapp.server.route('/download_sample/', methods=['GET'])
def download_sample():
    
    return send_from_directory(PathConfig.DATA_PATH, PathConfig.SAMPLE_DATA)