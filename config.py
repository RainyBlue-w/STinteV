import os

class PathConfig:

    ABS_ROOT_PATH = os.path.abspath(os.getcwd()) # 项目根目录绝对路径
    DATA_PATH = '/data1/share/omics-viewer/stintev'

class SecretConfig:
    SECRET_KEY = '4e17572a360111efaa763cecef387085'

class RouterConfig:

    # 合法pathname列表
    VALID_PATHNAME = [
        '/', 
        '/login',
        '/reset-password',
    ]
    
class NetConfig:
    
    SERVER_ADDRESS = 'http://10.86.60.31:8056'