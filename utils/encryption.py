import bcrypt

def pwd_encrypt(pwd: str):
    return bcrypt.hashpw(pwd.encode(encoding='utf-8'), bcrypt.gensalt())