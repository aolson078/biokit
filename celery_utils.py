from celery import Celery
import os


def make_celery(app):
    broker = app.config.get('CELERY_BROKER_URL', os.getenv('CELERY_BROKER_URL', 'redis://redis:6379/0'))
    backend = app.config.get('CELERY_RESULT_BACKEND', os.getenv('CELERY_RESULT_BACKEND', 'redis://redis:6379/0'))
    celery = Celery(app.import_name, broker=broker, backend=backend)
    celery.conf.update(app.config)

    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return super().__call__(*args, **kwargs)
    celery.Task = ContextTask
    return celery
