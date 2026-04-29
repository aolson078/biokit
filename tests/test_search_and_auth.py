import sys
sys.path.append('.')

import flask_bio_app as app_module
from flask_bio_app import create_app, db, Config
from models import User


class TestConfig(Config):
    TESTING = True
    SQLALCHEMY_DATABASE_URI = 'sqlite:///:memory:'
    WTF_CSRF_ENABLED = False
    CACHE_TYPE = 'SimpleCache'
    RATELIMIT_ENABLED = False


def test_search_normalizes_query_and_rejects_blank(monkeypatch):
    app = create_app(TestConfig)

    with app.app_context():
        db.create_all()
        user = User(username='searcher', email='search@example.com', raw_password='secret123')
        db.session.add(user)
        db.session.commit()
        user_id = user.id

    calls = []

    def fake_fetch_records(query):
        calls.append(query)
        return ['result-1']

    monkeypatch.setattr(app_module.utilities, 'fetch_records', fake_fetch_records)

    with app.test_client() as client:
        with client.session_transaction() as session:
            session['_user_id'] = str(user_id)
            session['_fresh'] = True

        blank = client.post('/search', json={'search': '   \n   '})
        assert blank.status_code == 400

        ok = client.post('/search', json={'search': '  BRCA1   human  '})
        assert ok.status_code == 200
        assert ok.get_json() == ['result-1']
        assert calls == ['BRCA1 human']
