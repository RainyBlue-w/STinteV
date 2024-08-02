from typing import Literal
from flask_login import UserMixin
from peewee import Model, SqliteDatabase, CharField
import os
from itsdangerous import URLSafeTimedSerializer
from stintev.config import PathConfig, SecretConfig
from stintev.utils import str2md5

db = SqliteDatabase(
    os.path.join(PathConfig.DATA_PATH, 'auth.db')
)

class User(UserMixin):
    pass

class UserAccount(Model):
    id = CharField(primary_key=True, unique=True)
    email = CharField(unique=True)
    password = CharField()
    user_level = CharField()
    
    class Meta:
        database = db
        table_name = 'user_account'

    @classmethod
    def add_one_user(cls, id, email, password, user_level: Literal['admin', 'normal']):

        return cls.insert_many([
            {
                'id': id,
                # password md5加密
                'password': str2md5(password),
                'email': email,
                'user_level': user_level
            }
        ]).on_conflict_replace().execute()
    
    @classmethod
    def query_user(cls, user_id):
        return list(UserAccount.select().where(UserAccount.id == user_id).dicts())
    
    # 获取密码重置token，10分钟有效
    def get_reset_token(self, expires_sec=600):
        s = URLSafeTimedSerializer(secret_key=SecretConfig.SECRET_KEY, expires_in=expires_sec)
        return s.dumps({'user_id': self.id}).decode('utf-8')
    
    # 验证密码重置token，返回UserAccount对象
    @staticmethod
    def verify_reset_token(token):
        s = URLSafeTimedSerializer(secret_key=SecretConfig.SECRET_KEY)
        try:
            user_id = s.loads(token)['user_id']
        except Exception as e:
            return None
        return UserAccount.query_user(user_id)[0]
    
db.create_tables([UserAccount])

if __name__ == '__main__':
    # 插入账号
    users = [
        {
            'id': 'admin',
            'username': 'admin',
            'email': '220222256@qq.com',
            'password': '72fad20d53965ba5ed36f3a280f7ff2d',
            'user_level': 'admin',
        },
        {
            'id': 'wuc',
            'username': 'wuc',
            'email': '1078497976@qq.com',
            'password': '72fad20d53965ba5ed36f3a280f7ff2d',
            'user_level': 'normal',
        }
    ]
    UserAccount.insert_many(users).execute()