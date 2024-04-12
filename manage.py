def deploy():
	from app import create_app, db
	from flask_migrate import upgrade, migrate, init, stamp

	app = create_app()
	# Tracks app level data, ensures functions have access to current app instance
	app.app_context().push()

	# migrate database to latest revision
	init()
	stamp()
	migrate()
	upgrade()

	db.create_all()


deploy()
