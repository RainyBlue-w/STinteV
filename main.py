# from fastapi import FastAPI
# from fastapi.middleware.wsgi import WSGIMiddleware
# from app import dash_app
# import argparse

# # args

# parser = argparse.ArgumentParser(description='Run the 3Dviewer app')
# parser.add_argument('--dir-data', type=str, default='data', help='Directory containing the h5ad files'),
# parser.add_argument('--palette', type=str, default='data/celltype_cmap.csv', help='Celltype color palette to use'),
# parser.add_argument('--dir-cache', type=str, default='cache', help='Directory to store the cache files'),
# parser.add_argument('--port', type=int, default=8000, help='Port to run the app on')
# parser.add_argument('--host', type=str, default='0.0.0.0', help='Host to run the app on')
# parser.add_argument('--path-mount', type=str, default='/', help='Path to mount the app on'),
# args = parser.parse_args()

# # app
# app = FastAPI()
# app.mount(args.path_mount, WSGIMiddleware(dash_app.server))

# if __name__ == "__main__":
#     import uvicorn
#     uvicorn.run(
#         app, 
#         port=args.port,
#         host=args.host,
#         log_level="info"
#     )
