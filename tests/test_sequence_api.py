import re

from flask_bio_app import Config, create_app, db
from models import Record, User


class SequenceApiConfig(Config):
    TESTING = True
    SQLALCHEMY_DATABASE_URI = "sqlite:///:memory:"
    SESSION_COOKIE_SECURE = False
    CACHE_TYPE = "SimpleCache"
    RATELIMIT_ENABLED = False
    SEQUENCE_WORKSPACE_ENABLED = True
    SEQUENCE_WORKSPACE_CTA_ENABLED = False
    VITE_DEV_SERVER_URL = "http://localhost:5173"
    CELERY_TASK_ALWAYS_EAGER = True
    CELERY_TASK_EAGER_PROPAGATES = True


def create_user(app, username="sequence-user", email="sequence@example.com"):
    with app.app_context():
        db.create_all()
        user = User(username=username, email=email, raw_password="secret123")
        db.session.add(user)
        db.session.commit()
        return user.id


def authenticate(client, user_id):
    with client.session_transaction() as session:
        session["_user_id"] = str(user_id)
        session["_fresh"] = True


def csrf_token(client):
    page = client.get("/sequence-wizard")
    assert page.status_code == 200
    match = re.search(rb'<meta name="csrf-token" content="([^"]+)"', page.data)
    assert match
    return match.group(1).decode()


def test_api_authentication_and_csrf_are_json_contracts():
    app = create_app(SequenceApiConfig)
    user_id = create_user(app)
    with app.test_client() as client:
        unauthenticated = client.get("/api/v1/analyses")
        assert unauthenticated.status_code == 401
        assert unauthenticated.get_json()["error"]["code"] == "authentication_required"

        authenticate(client, user_id)
        missing_csrf = client.post("/api/v1/nucleotides/import", json={"text": "ATGTAA", "molecule_type": "dna"})
        assert missing_csrf.status_code == 400
        assert missing_csrf.get_json()["error"]["code"] == "csrf_failed"

        token = csrf_token(client)
        imported = client.post(
            "/api/v1/nucleotides/import",
            json={"text": "ATGTAA", "molecule_type": "dna"},
            headers={"X-CSRFToken": token},
        )
        assert imported.status_code == 200
        assert imported.get_json()["records"][0]["sequence"] == "ATGTAA"


def test_inline_summary_does_not_write_and_save_is_explicit_and_owner_scoped():
    app = create_app(SequenceApiConfig)
    owner_id = create_user(app)
    with app.app_context():
        other = User(username="other-user", email="other@example.com", raw_password="secret123")
        db.session.add(other)
        db.session.commit()
        other_id = other.id

    with app.test_client() as client:
        authenticate(client, owner_id)
        token = csrf_token(client)
        summary = client.post(
            "/api/v1/analyses/summary",
            json={"record": {"source": "manual", "sequence": "ATGTAA", "molecule_type": "dna"}},
            headers={"X-CSRFToken": token},
        )
        assert summary.status_code == 200
        with app.app_context():
            assert Record.query.count() == 0

        saved = client.post(
            "/api/v1/records",
            json={"source": "manual", "sequence": "ATGTAA", "molecule_type": "dna", "user_label": "Example"},
            headers={"X-CSRFToken": token},
        )
        assert saved.status_code == 201
        record_id = saved.get_json()["record"]["id"]

    with app.test_client() as other_client:
        authenticate(other_client, other_id)
        hidden = other_client.get(f"/api/v1/records/{record_id}")
        assert hidden.status_code == 404


def test_legacy_process_sequence_is_authenticated_gone_without_writing():
    app = create_app(SequenceApiConfig)
    user_id = create_user(app)
    with app.test_client() as client:
        authenticate(client, user_id)
        token = csrf_token(client)
        response = client.post("/process_sequence", json={"sequence": "ATGTAA"}, headers={"X-CSRFToken": token})
        assert response.status_code == 410
        assert response.get_json()["error"]["code"] == "legacy_save_retired"
        with app.app_context():
            assert Record.query.count() == 0
