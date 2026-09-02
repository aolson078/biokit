"""Schema-bearing baseline for the shipped BioKit SQLite database.

Revision ID: 0001_shipped_schema_baseline
Revises:
"""

from alembic import op
import sqlalchemy as sa

revision = "0001_shipped_schema_baseline"
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "user",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("username", sa.String(100), nullable=False, unique=True),
        sa.Column("email", sa.String(120), nullable=False, unique=True),
        sa.Column("password", sa.String(300), nullable=False),
        sa.Column("role", sa.String(40)),
        sa.Column("view_reports", sa.Boolean()),
        sa.Column("delete_reports", sa.Boolean()),
        sa.Column("print_reports", sa.Boolean()),
        sa.Column("change_reports", sa.Boolean()),
    )
    op.create_table(
        "report",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("employee_id", sa.Integer(), sa.ForeignKey("user.id")),
        sa.Column("nucleotide_ids", sa.JSON()),
        sa.Column("organisms", sa.JSON()),
        sa.Column("nucleotides", sa.JSON()),
        sa.Column("phylo_tree", sa.String(50)),
        sa.Column("dot_line_graph", sa.JSON()),
        sa.Column("heat_map", sa.JSON()),
        sa.Column("bar_chart", sa.String(50)),
    )
    op.create_table(
        "record",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("nucleotide_id", sa.String(20), unique=True),
        sa.Column("organism", sa.String(80)),
        sa.Column("gene_info", sa.String(100)),
        sa.Column("nucleotides", sa.Text()),
        sa.Column("gc_content", sa.Float()),
        sa.Column("amino_acids", sa.Text()),
        sa.Column("siRNA", sa.Text()),
        sa.Column("sense_similarity", sa.Float()),
        sa.Column("melting_temp", sa.Float()),
        sa.Column("molecular_weight", sa.Float()),
        sa.Column("hydrophobicity", sa.Float()),
        sa.Column("secondary_structure_prediction", sa.String(100)),
        sa.Column("employee_id", sa.Integer(), sa.ForeignKey("user.id"), nullable=False),
        sa.Column("report_id", sa.Integer(), sa.ForeignKey("report.id")),
    )
    op.create_table(
        "report_record",
        sa.Column("report_id", sa.Integer(), sa.ForeignKey("report.id"), primary_key=True),
        sa.Column("record_id", sa.Integer(), sa.ForeignKey("record.id"), primary_key=True),
    )


def downgrade():
    raise RuntimeError("Baseline downgrade is destructive; restore the verified backup instead.")
