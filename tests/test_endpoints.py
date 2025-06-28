import sys
sys.path.append('.')
from flask_bio_app import create_app


def test_index_route():
    app = create_app()
    with app.test_client() as client:
        resp = client.get('/')
        assert resp.status_code == 200



def test_protected_reports_requires_login():
    app = create_app()
    with app.test_client() as client:
        resp = client.get('/api/v1/reports')
        assert resp.status_code in (302, 401)
