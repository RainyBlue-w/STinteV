import hashlib

def str2md5(string):
    return hashlib.md5(string.encode('utf-8')).hexdigest()