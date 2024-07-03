import os

class PathConfig:

    ABS_ROOT_PATH = os.path.abspath(os.getcwd()) # 项目根目录绝对路径
    DATA_PATH = '/rad/share/omics-viewer/stintev'


class RouterConfig:

    # 合法pathname列表
    VALID_PATHNAME = [
        '/', '/login'
    ]