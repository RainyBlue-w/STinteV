from typing import Literal, List
from flask_login import UserMixin
from peewee import Model, SqliteDatabase, CharField
import os
from itsdangerous import URLSafeTimedSerializer
from stintev.config import PathConfig, SecretConfig
from stintev.utils import str2md5

# database for user accounts
db = SqliteDatabase(
    os.path.join(PathConfig.DATA_PATH, 'auth.db')
)

class User(Model, UserMixin):
    
    id = CharField(primary_key=True, unique=True)
    email = CharField(unique=True)
    password = CharField()
    user_level = CharField()
    
    class Meta:
        database = db
        table_name = 'user_account'
    
    # 获取密码重置token，10分钟有效
    def get_reset_token(self):
        s = URLSafeTimedSerializer(secret_key=SecretConfig.SECRET_KEY)
        return s.dumps({'user_id': self.id})
    
    # 验证密码重置token，返回User对象
    @staticmethod
    def verify_reset_token(token):
        s = URLSafeTimedSerializer(secret_key=SecretConfig.SECRET_KEY)
        try:
            user_id = s.loads(token, max_age=600)['user_id']
        except Exception as e:
            return None
        return User.query_by_id(user_id)[0]
    
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
    def query_by_id(cls, user_id):
        return list(User.select().where(User.id == user_id))
    
    @classmethod
    def query_by_email(cls, email):
        return list(User.select().where(User.email == email))
    
    @classmethod
    def update_password(cls, user_id, new_password):
        try:
            User.query_by_id(user_id)
            res=(User
                 .update({User.password: str2md5(new_password)})
                 .where(User.id == user_id)
                 .execute())
            return True
        except Exception as e:
            return False

db.create_tables([User])

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
    User.insert_many(users).execute()
    User