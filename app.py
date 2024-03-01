from flask import Flask

# Application generator ------------------------------------------------------------------------------------------------
def create_app():
	app = Flask(__name__)
	app.secret_key = 'your-secret-key'
	return app


if __name__ == '__main__':
	app = create_app()
	app.run()
