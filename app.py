from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_bcrypt import Bcrypt
from flask_migrate import Migrate
from flask_login import LoginManager

# Session manager ------------------------------------------------------------------------------------------------------
# session management, handles login and current user information and permissions
login_manager = LoginManager()
# strong = Session will be invalid if the users IP address changes during session
login_manager.session_protection = "strong"
# sends user to login if attempting to access resources while not authenticated
login_manager.login_view = "login"
login_manager.login_message_category = "info"




# Database -------------------------------------------------------------------------------------------------------------
# handles SQL queries and database interactions
db = SQLAlchemy()
# handles updates in db schema
migrate = Migrate()
# handles encryption for passwords
bcrypt = Bcrypt()


# Application generator ------------------------------------------------------------------------------------------------
def create_app():
	# creates instance of Flask application with current files directory as root
	app = Flask(__name__)
	# set key to sign session cookies (in production should save as env variable on host machine or other means of secrecy)
	app.secret_key = "was-software"
	# path where Flask looks to connect to database
	app.config['SQLALCHEMY_DATABASE_URI'] = "sqlite:///database.db"
	# supress tracking changes in db, adds unnecessary overhead
	app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

	# initiate instances in application context
	login_manager.init_app(app)
	db.init_app(app)
	migrate.init_app(app, db)
	bcrypt.init_app(app)

	#removed
	#return app

	const = "THIS IS FOREVER"
	
	#Will
	#added
	with app.app_context():
		db.create_all()
	return app

#removed
#if __name__ == '__main__':
	#app = create_app()
	#app.run()
