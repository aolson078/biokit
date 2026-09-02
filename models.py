"""Database models for the BioKit application."""

from datetime import datetime

from flask_bcrypt import generate_password_hash, check_password_hash
from flask_login import UserMixin
from sqlalchemy import text
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

    # Keep the safer Python attribute while matching the shipped database column.
    password_hash = db.Column("password", db.String(300), nullable=False)

    role            = db.Column(db.String(40), default="employee")
    view_reports    = db.Column(db.Boolean, default=False)
    delete_reports  = db.Column(db.Boolean, default=False)
    print_reports   = db.Column(db.Boolean, default=False)
    change_reports  = db.Column(db.Boolean, default=False)

    # ── convenience setters/getters ──────────────────────────────
    def set_password(self, raw_password: str) -> None:
        self.password_hash = generate_password_hash(raw_password).decode("utf-8")

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

    __table_args__ = (
        db.CheckConstraint("source IN ('ncbi','manual','legacy')", name="record_source"),
        db.CheckConstraint("molecule_type IN ('dna','rna','unknown')", name="record_molecule_type"),
        db.CheckConstraint("sequence_alphabet IN ('dna','rna','neutral','unknown')", name="record_sequence_alphabet"),
        db.Index(
            "uq_record_employee_source_accession",
            "employee_id",
            "source_accession",
            unique=True,
            sqlite_where=text("source_accession IS NOT NULL"),
        ),
    )

    id = db.Column(db.Integer, primary_key=True)
    nucleotide_id = db.Column(db.String(64))
    source = db.Column(db.String(16), nullable=False, default="legacy")
    source_accession = db.Column(db.String(64))
    source_title = db.Column(db.Text)
    user_label = db.Column(db.String(120))
    source_retrieved_at = db.Column(db.DateTime(timezone=True))
    source_updated_at = db.Column(db.DateTime(timezone=True))
    molecule_type = db.Column(db.String(8), nullable=False, default="unknown")
    sequence_alphabet = db.Column(db.String(8), nullable=False, default="unknown")
    sequence_length = db.Column(db.Integer)
    created_at = db.Column(db.DateTime(timezone=True), nullable=False, default=datetime.utcnow)
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
    entropy_line_graphs = db.Column(db.JSON)
    kmer_histograms = db.Column(db.JSON)
    cumulative_gc_skew_graphs = db.Column(db.JSON)
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
