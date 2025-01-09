import flask
from werkzeug.serving import run_simple


def create_flask_app():
    app = flask.Flask(__name__)

    @app.route("/")
    def home():
        return "Hello, Flask!"

    return app