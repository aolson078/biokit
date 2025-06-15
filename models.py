from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin
from app import db

class User(UserMixin, db.Model):
    __tablename__ = "user"

    id       = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(100), unique=True, nullable=False)
    email    = db.Column(db.String(120), unique=True, nullable=False)

    password_hash = db.Column(db.String(256), nullable=False)

    role            = db.Column(db.String(40), default="employee")
    view_reports    = db.Column(db.Boolean, default=False)
    delete_reports  = db.Column(db.Boolean, default=False)
    print_reports   = db.Column(db.Boolean, default=False)
    change_reports  = db.Column(db.Boolean, default=False)

    # ── convenience setters/getters ──────────────────────────────
    def set_password(self, raw_password: str) -> None:
        self.password_hash = generate_password_hash(raw_password)

    def check_password(self, raw_password: str) -> bool:
        return check_password_hash(self.password_hash, raw_password)

    # optional: accept raw_password in __init__ for one-liners
    def __init__(self, username: str, email: str, raw_password: str, role: str = "employee"):
        self.username = username
        self.email    = email
        self.set_password(raw_password)   # hashes internally
        self.role     = role

    def __repr__(self):
        return f"<User {self.username}>"
