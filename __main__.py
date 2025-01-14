import argparse
from colorama import init
from stintev.app import dashapp
from stintev.config import NetConfig

import uvicorn
from fastapi import FastAPI
from fastapi.middleware.wsgi import WSGIMiddleware

parser = argparse.ArgumentParser(description='Run the 3Dviewer app')
parser.add_argument('--host', default='0.0.0.0', required=False, help='Host to run the app on')
parser.add_argument('--port', default=443, required=False, help='Port to run the app on')
parser.add_argument('--init-db', action='store_true', help='Initialize the database')
parser.add_argument('--dev', action='store_true', help='Run the app in dev mode')
parser.add_argument('--debug', action='store_true', help='Run the app in debug mode')
args, unknown = parser.parse_known_args()

if args.init_db:
    from stintev.models.init_db import init_db
    init_db()
    
elif args.dev:
    dashapp.register_celery_tasks()
    dashapp.run(
        host = args.host,
        port = int(args.port),
        debug = args.debug,
        use_reloader = False,
    )
else:
    app = FastAPI()
    app.mount("/", WSGIMiddleware(dashapp.server))
        
if __name__ == "__main__":
    uvicorn.run(
        app,
        host = args.host,
        port = int(args.port),
        ssl_keyfile = NetConfig.SSL_KEYFILE,
        ssl_certfile = NetConfig.SSL_CRTFILE
    )