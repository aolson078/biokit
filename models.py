"""Database models for the BioKit application."""

from datetime import datetime

from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin
from flask_bio_app import db

# Association table between reports and records
report_record = db.Table(
    "report_record",
    db.Column("report_id", db.Integer, db.ForeignKey("report.id"), primary_key=True),
    db.Column("record_id", db.Integer, db.ForeignKey("record.id"), primary_key=True),
)

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

    # ----- permission helpers -------------------------------------------------
    def can_view_report(self, report) -> bool:
        return self.view_reports or report.employee_id == self.id

    def can_delete_report(self, report) -> bool:
        return self.delete_reports or report.employee_id == self.id

    def can_print_report(self, report) -> bool:
        return self.print_reports or report.employee_id == self.id

    def can_change_report(self, report) -> bool:
        return self.change_reports or report.employee_id == self.id


class Record(db.Model):
    """Individual sequence analysis results."""

    id = db.Column(db.Integer, primary_key=True)
    nucleotide_id = db.Column(db.String(20), unique=True)
    organism = db.Column(db.String(80))
    gene_info = db.Column(db.String(100))
    nucleotides = db.Column(db.Text)
    gc_content = db.Column(db.Float)
    amino_acids = db.Column(db.Text)
    siRNA = db.Column(db.Text)
    sense_similarity = db.Column(db.Float)
    melting_temp = db.Column(db.Float)
    molecular_weight = db.Column(db.Float)
    hydrophobicity = db.Column(db.Float)
    secondary_structure_prediction = db.Column(db.String(100))

    employee_id = db.Column(db.Integer, db.ForeignKey("user.id"), nullable=False)

    def __repr__(self):
        return f"<Record {self.nucleotide_id}>"


class Report(db.Model):
    """Compiled report containing multiple records."""

    id = db.Column(db.Integer, primary_key=True)
    employee_id = db.Column(db.Integer, db.ForeignKey("user.id"))

    nucleotide_ids = db.Column(db.JSON)
    organisms = db.Column(db.JSON)
    nucleotides = db.Column(db.JSON)

    phylo_tree = db.Column(db.String(50))
    dot_line_graph = db.Column(db.JSON)
    heat_map = db.Column(db.JSON)
    bar_chart = db.Column(db.String(50))
    gc_line_graphs = db.Column(db.JSON)
    gc_skew_graphs = db.Column(db.JSON)
    nuc_pie_charts = db.Column(db.JSON)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    employee = db.relationship("User", backref="reports")
    associated_records = db.relationship(
        "Record",
        secondary=report_record,
        lazy="dynamic",
        backref=db.backref("reports", lazy="dynamic"),
    )

    def __repr__(self):
        return f"<Report {self.id}>"
