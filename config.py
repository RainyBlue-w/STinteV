import os

class PathConfig:

    ABS_ROOT_PATH = os.path.abspath(os.getcwd()) # 项目根目录绝对路径
    DATA_PATH = '/data/data_stintev/'

class SecretConfig:
    SECRET_KEY = '4e17572a360111efaa763cecef387085'

class RouterConfig:

    # 合法pathname列表
    VALID_PATHNAME = [
        '/', 
        # '/login',
        # '/reset-password',
    ]
    
class NetConfig:
    
    SERVER_ADDRESS = 'https://stintev.site'
    SSL_KEYFILE = '/home/wuc/stintev.site_other/stintev.site.key'
    SSL_CRTFILE = '/home/wuc/stintev.site_other/stintev.site_bundle.pem'