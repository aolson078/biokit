from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_bcrypt import Bcrypt
from flask_migrate import Migrate

from flask_login import (
	UserMixin,
	login_user,
	LoginManager,
	current_user,
	logout_user,
	login_required,
)

# Session manager ------------------------------------------------------------------------------------------------------
# Login_manager object handles auth and session management
login_manager = LoginManager()
# Strong = Session will be invalid if the users IP address changes during session
login_manager.session_protection = "strong"
# Sends user to login if attempting to access resources while not authenticated
login_manager.login_view = "login"
login_manager.login_message_category = "info"

# Database -------------------------------------------------------------------------------------------------------------
db = SQLAlchemy()
migrate = Migrate()
bcrypt = Bcrypt()


# Application generator ------------------------------------------------------------------------------------------------
def create_app():
	app = Flask(__name__)
	app.secret_key = "was-software"
	app.config['SQLALCHEMY_DATABASE_URI'] = "sqlite:///database.db"
	app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = True

	login_manager.init_app(app)
	db.init_app(app)
	migrate.init_app(app, db)
	bcrypt.init_app(app)
	print("Ran create app successfully")

	return app


if __name__ == '__main__':
	app = create_app()
	app.run()
