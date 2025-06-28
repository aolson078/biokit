import sys
sys.path.append('.')
from flask_bio_app import create_app

def test_create_app():
    app = create_app()
    assert app is not None
    assert hasattr(app, 'celery')
