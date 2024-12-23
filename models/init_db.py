from flask_login import UserMixin
from peewee import Model, SqliteDatabase, CharField
import os
from itsdangerous import URLSafeTimedSerializer
from stintev.config import PathConfig

# database for user accounts
db = SqliteDatabase(
    os.path.join(PathConfig.DATA_PATH, 'auth.db')
)

class User(Model, UserMixin):
    
    id = CharField(primary_key=True, unique=True)
    
    class Meta:
        database = db
        table_name = 'user_account'
    
    @classmethod
    def add_one_user(cls, id):

        return cls.insert_many([
            {
                'id': id,
            }
        ]).on_conflict_replace().execute()
    
    @classmethod
    def query_by_id(cls, user_id):
        return list(User.select().where(User.id == user_id))



def init_db():    
    db.create_tables([User])