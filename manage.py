def deploy():
	from app import create_app, db
	from flask_migrate import upgrade, migrate, init, stamp

	# create instance of flask app
	app = create_app()

	# tracks app level data, ensures functions have access to current app instance
	app.app_context().push()

	# migrate database to latest revision and initialize a new migration directory
	init()

	# stamp the migration directory with the current revision
	stamp()

	# autogenerate a new revision file
	migrate()

	# upgrade to the latest revision
	upgrade()

	# create any tables defined in the models that don't exist yet
	db.create_all()


deploy()
