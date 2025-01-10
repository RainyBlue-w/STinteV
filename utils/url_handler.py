import urllib.parse

def decompress_list(param: str):
    '''
    param: 'param1%2C%param2%2C%param3'  # %2C% -> ','
    '''

    list = param.split('%2C%')
    
    return list
    
    
    