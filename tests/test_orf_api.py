import sys
sys.path.append('.')

from flask_bio_app import create_app, db, Config
from models import User


class TestConfig(Config):
    TESTING = True
    SQLALCHEMY_DATABASE_URI = 'sqlite:///:memory:'
    WTF_CSRF_ENABLED = False
    CACHE_TYPE = 'SimpleCache'
    RATELIMIT_ENABLED = False


def _login_session(client, user_id):
    with client.session_transaction() as session:
        session['_user_id'] = str(user_id)
        session['_fresh'] = True


def test_orf_api_accepts_fasta_and_returns_metadata():
    app = create_app(TestConfig)
    with app.app_context():
        db.create_all()
        user = User(username='orfuser', email='orf@example.com', raw_password='secret123')
        db.session.add(user)
        db.session.commit()
        user_id = user.id

    with app.test_client() as client:
        _login_session(client, user_id)
        payload = {'sequence': '>seq1\nCCCATGAAAAAAAAAAAAAAAAAAAAAAAAAAATAAGGG', 'min_length': 30}
        response = client.post('/api/orf-finder', json=payload)

        assert response.status_code == 200
        data = response.get_json()
        assert data['sequence_length'] == 39
        assert data['count'] >= 1
        first = data['orfs'][0]
        assert {'start', 'end', 'length_nt', 'length_aa', 'protein', 'strand', 'frame'} <= set(first)
        assert first['start'] >= 1


def test_orf_api_rejects_invalid_sequence():
    app = create_app(TestConfig)
    with app.app_context():
        db.create_all()
        user = User(username='orfbad', email='orfbad@example.com', raw_password='secret123')
        db.session.add(user)
        db.session.commit()
        user_id = user.id

    with app.test_client() as client:
        _login_session(client, user_id)
        response = client.post('/api/orf-finder', json={'sequence': 'ATGBZX', 'min_length': 30})
        assert response.status_code == 400
        assert 'Invalid DNA sequence' in response.get_json()['error']
