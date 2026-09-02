"""Add sequence provenance and user-scoped accession identity.

Revision ID: 0002_sequence_workspace_metadata
Revises: 0001_shipped_schema_baseline
"""

from alembic import op
import sqlalchemy as sa

revision = "0002_sequence_workspace_metadata"
down_revision = "0001_shipped_schema_baseline"
branch_labels = None
depends_on = None

NAMING_CONVENTION = {
    "ix": "ix_%(column_0_label)s",
    "uq": "uq_%(table_name)s_%(column_0_name)s",
    "ck": "ck_%(table_name)s_%(constraint_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
    "pk": "pk_%(table_name)s",
}


def upgrade():
    with op.batch_alter_table("record", recreate="always", naming_convention=NAMING_CONVENTION) as batch:
        batch.alter_column("nucleotide_id", existing_type=sa.String(20), type_=sa.String(64))
        batch.drop_constraint("uq_record_nucleotide_id", type_="unique")
        batch.add_column(sa.Column("source", sa.String(16), nullable=False, server_default="legacy"))
        batch.add_column(sa.Column("source_accession", sa.String(64)))
        batch.add_column(sa.Column("source_title", sa.Text()))
        batch.add_column(sa.Column("user_label", sa.String(120)))
        batch.add_column(sa.Column("source_retrieved_at", sa.DateTime(timezone=True)))
        batch.add_column(sa.Column("source_updated_at", sa.DateTime(timezone=True)))
        batch.add_column(sa.Column("molecule_type", sa.String(8), nullable=False, server_default="unknown"))
        batch.add_column(sa.Column("sequence_alphabet", sa.String(8), nullable=False, server_default="unknown"))
        batch.add_column(sa.Column("sequence_length", sa.Integer()))
        batch.add_column(sa.Column("created_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.text("CURRENT_TIMESTAMP")))
        batch.create_check_constraint("record_source", "source IN ('ncbi','manual','legacy')")
        batch.create_check_constraint("record_molecule_type", "molecule_type IN ('dna','rna','unknown')")
        batch.create_check_constraint("record_sequence_alphabet", "sequence_alphabet IN ('dna','rna','neutral','unknown')")

    op.execute("UPDATE record SET sequence_length = length(nucleotides) WHERE nucleotides IS NOT NULL")
    op.create_index(
        "uq_record_employee_source_accession",
        "record",
        ["employee_id", "source_accession"],
        unique=True,
        sqlite_where=sa.text("source_accession IS NOT NULL"),
    )


def downgrade():
    raise RuntimeError("Use the verified pre-migration backup; global uniqueness may no longer be restorable.")
