"""Canonical Celery application and registered BioKit tasks."""

from flask_bio_app import create_app

app = create_app()
celery = app.celery
