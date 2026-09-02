import sys
sys.path.append('.')

from flask_bio_app import create_app, db, Config
from models import User
from sequence_workspace.ncbi import NcbiClient


class TestConfig(Config):
    TESTING = True
    SQLALCHEMY_DATABASE_URI = 'sqlite:///:memory:'
    SESSION_COOKIE_SECURE = False
    WTF_CSRF_ENABLED = False
    CACHE_TYPE = 'SimpleCache'
    RATELIMIT_ENABLED = False
    NCBI_EMAIL = 'test@example.com'


def test_user_model_maps_hash_to_shipped_password_column():
    assert User.password_hash.property.columns[0].name == 'password'


def test_login_accepts_valid_credentials():
    app = create_app(TestConfig)

    with app.app_context():
        db.create_all()
        db.session.add(User(username='loginuser', email='login@example.com', raw_password='secret123'))
        db.session.commit()

    with app.test_client() as client:
        response = client.post('/login', data={
            'email': 'login@example.com',
            'password': 'secret123',
        })

        with client.session_transaction() as session:
            assert session['_user_id']

        authenticated_page = client.get(response.headers['Location'])

    assert response.status_code == 302
    assert response.headers['Location'].endswith('/')
    assert authenticated_page.status_code == 200


def test_login_rejects_invalid_password():
    app = create_app(TestConfig)

    with app.app_context():
        db.create_all()
        db.session.add(User(username='loginuser', email='login@example.com', raw_password='secret123'))
        db.session.commit()

    with app.test_client() as client:
        response = client.post('/login', data={
            'email': 'login@example.com',
            'password': 'wrong-password',
        })

        with client.session_transaction() as session:
            assert '_user_id' not in session

    assert response.status_code == 200
    assert b'Invalid email or password' in response.data


def test_search_normalizes_query_and_rejects_blank(monkeypatch):
    app = create_app(TestConfig)

    with app.app_context():
        db.create_all()
        user = User(username='searcher', email='search@example.com', raw_password='secret123')
        db.session.add(user)
        db.session.commit()
        user_id = user.id

    calls = []

    def fake_search(self, query, organism=None, page=1, page_size=10):
        calls.append(query)
        return {'results': [{
            'accession_version': 'NM_1.1',
            'title': 'result-1',
            'organism': 'Testus',
            'length': 42,
        }]}

    monkeypatch.setattr(NcbiClient, 'search', fake_search)

    with app.test_client() as client:
        with client.session_transaction() as session:
            session['_user_id'] = str(user_id)
            session['_fresh'] = True

        blank = client.post('/search', json={'search': '   \n   '})
        assert blank.status_code == 400

        ok = client.post('/search', json={'search': '  BRCA1   human  '})
        assert ok.status_code == 200
        assert ok.get_json() == [{
            'title': 'result-1',
            'description': 'NM_1.1',
            'organism': 'Testus',
            'length': 42,
        }]
        assert calls == ['BRCA1 human']
