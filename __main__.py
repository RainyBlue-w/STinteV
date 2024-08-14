import argparse
from stintev.app import dashapp

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run the 3Dviewer app')
    parser.add_argument('--host', default='0.0.0.0', required=False, help='Host to run the app on')
    parser.add_argument('--port', default='8055', required=False, help='Port to run the app on')
    parser.add_argument('--debug', action='store_true', required=False, help='Run app on debug mode')
    args = parser.parse_args()
    
    dashapp.run(
        host=args.host,
        port=args.port,
        debug=args.debug
    )