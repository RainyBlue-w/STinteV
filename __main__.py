import argparse
from colorama import init
from stintev.app import dashapp
from stintev.config import NetConfig

import uvicorn
from fastapi import FastAPI
from fastapi.middleware.wsgi import WSGIMiddleware

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run the 3Dviewer app')
    parser.add_argument('--host', default='0.0.0.0', required=False, help='Host to run the app on')
    parser.add_argument('--port', default=443, required=False, help='Port to run the app on')
    parser.add_argument('--init-db', action='store_true', help='Initialize the database')
    parser.add_argument('--dev', action='store_true', help='Run the app in dev mode')
    parser.add_argument('--debug', action='store_true', help='Run the app in debug mode')
    args = parser.parse_args()
    
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
        dashapp.index_string = '''
        <!DOCTYPE html>
        <html>
            <script defer src="http://localhost:3000/random-string.js" data-website-id="1d5775bf-1485-4b2e-949e-1e23c8d9625d"></script>
            <head>
                {%metas%}
                <title>{%title%}</title>
                {%favicon%}
                {%css%}
            </head>
            <body>
                {%app_entry%}
                <footer>
                    {%config%}
                    {%scripts%}
                    {%renderer%}
                </footer>
                <script type="text/javascript" id="cookiebanner"
                    src="https://cdn.jsdelivr.net/gh/dobarkod/cookie-banner@1.2.2/dist/cookiebanner.min.js"
                    data-height="36px" data-position="bottom" data-font-size="18px"
                    data-close-text="Got it!"
                    data-message="We use cookies to preserve your session state, by using this website you agree to our use of cookies.">
                </script>
            </body>
        </html>
        '''
        app = FastAPI()
        app.mount("/", WSGIMiddleware(dashapp.server))
        uvicorn.run(
            app,
            host = args.host,
            port = int(args.port),
            ssl_keyfile = NetConfig.SSL_KEYFILE,
            ssl_certfile = NetConfig.SSL_CRTFILE
        )