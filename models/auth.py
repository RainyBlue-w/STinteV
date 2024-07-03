from typing import Literal
from flask_login import UserMixin
from peewee import Model, SqliteDatabase, CharField
import os
from stintev.config import PathConfig
from stintev.utils import str2md5

db = SqliteDatabase(
    os.path.join(PathConfig.DATA_PATH, 'auth.db')
)

class User(UserMixin):
    pass

class UserAccount(Model):
    id = CharField(primary_key=True, unique=True)
    username = CharField(unique=True)
    email = CharField()
    password = CharField()
    user_level = CharField()
    
    class Meta:
        database = db
        table_name = 'user_account'

    @classmethod
    def add_one_user(cls, username, email, password, user_level: Literal['admin', 'normal']):

        return cls.insert_many([
            {
                'id': username,
                'username': username,
                # password md5加密
                'password': str2md5(password),
                'email': email,
                'user_level': user_level
            }
        ]).on_conflict_replace().execute()
        
    @classmethod
    def query_user(cls, username):
        return list(UserAccount.select().where(UserAccount.username== username).dicts())
        
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