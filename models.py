from functools import wraps

from Bio import Entrez
from flask import redirect, url_for

from custom_exceptions import InvalidQueryException
from app import create_app, db
from flask_login import UserMixin, current_user


# UserMixin gives us premade methods for authentication and user management
class User(UserMixin, db.Model):
	__tablename__ = 'user'

	id = db.Column(db.Integer, primary_key=True)
	username = db.Column(db.String(100), unique=True, nullable=False)
	email = db.Column(db.String(120), unique=True, nullable=False)
	password = db.Column(db.String(300), nullable=False, unique=True)
	role = db.Column(db.String(40), nullable=True)

	# looks up user in db and adds document to list
	def __init__(self, username, email, password, role='employee'):
		self.report_list = None
		self.username = username
		self.email = email
		self.password = password
		self.role = role

	def add_report(self, report):
		if self.report_list is None:
			self.report_list = []
		self.report_list.append(report)
		db.session.commit()

	def __repr__(self):
		return '<User: %r>' % self.username


class Record(db.Model):
	__tablename__ = 'record'

	id = db.Column(db.Integer, primary_key=True)
	nucleotide_id = db.Column(db.String(20), unique=True)
	organism = db.Column(db.String(80))
	gene_info = db.Column(db.String(100))
	nucleotides = db.Column(db.Text)

	# single record calculations
	gc_content = db.Column(db.Float)
	amino_acids = db.Column(db.Text)
	siRNA = db.Column(db.Text)
	sense_similarity = db.Column(db.Float)
	melting_temp = db.Column(db.Float)
	molecular_weight = db.Column(db.Float)
	hydrophobicity = db.Column(db.Float)
	secondary_structure_prediction = db.Column(db.String(100))

	report_id = db.Column(db.Integer, db.ForeignKey('report.id'))
	reports = db.relationship('Report', secondary='report_record', backref='associated_records')


report_record = db.Table('report_record',
                         db.Column('report_id', db.Integer, db.ForeignKey('report.id'), primary_key=True),
                         db.Column('record_id', db.Integer, db.ForeignKey('record.id'), primary_key=True)
                         )


# the report class represents the final product. It will contain the computed data from the bio processes
class Report(db.Model):
	__tablename__ = 'report'
	id = db.Column(db.Integer, primary_key=True)
	employee_id = db.Column(db.Integer, nullable=False)
	# holds ids of all nuc strings used in calculations
	nucleotide_ids = db.Column(db.JSON)
	# holds name of each organism in report
	organisms = db.Column(db.JSON)
	nucleotides = db.Column(db.JSON)
	phylo_tree = db.Column(db.String(50), nullable=True)
	dot_line_graph = db.Column(db.JSON)
	heat_map = db.Column(db.JSON)
	bar_chart = db.Column(db.String(50), nullable=True)
	records = db.relationship('Record', secondary='report_record', backref='associated_reports')


def is_manager(func):
	@wraps(func)
	def authenticate_manager_view(*args, **kwargs):
		if not current_user.is_authenticated or current_user.role != 'manager':
			return redirect(url_for('denied'))
		return func(*args, **kwargs)
	return authenticate_manager_view


def is_admin(func):
	@wraps(func)
	def authenticate_admin_view(*args, **kwargs):
		if not current_user.is_authenticated or current_user.role != 'admin':
			return redirect(url_for('denied'))
		return func(*args, **kwargs)
	return authenticate_admin_view


# queries the selected database for term and returns the record with the nucleotide string
def fetch_records(query):
	Entrez.email = "aolson078@gmail.com"
	Entrez.apikey = "4a1d5a80996f3691a67335ded0b79299c708"
	# nucleotide, gene, or protein for db
	database = "nucleotide"
	# number of desired responses per query
	return_max = 5
	# fetch nucleotide record related to term
	IDs = Entrez.read(Entrez.esearch(db=database, term=query, field="Organism", retmax=return_max))["IdList"]
	records = []
	for ID in IDs:
		handle = Entrez.efetch(db=database, id=ID, rettype="fasta", retmod="text")
		fasta_record = handle.read()

		# Parse the FASTA record to extract the title and description
		lines = fasta_record.split("\n")
		title = lines[0].strip(">")  # Extract the title from the first line
		description = "\n".join(lines[1:])  # Join the remaining lines as the description

		record = {
			"id": ID,
			"title": title,
			"description": description
		}
		records.append(record)

	if len(records) == 0:
		raise InvalidQueryException(query, database)

	return records


if __name__ == '__main__':
	app = create_app()
	with app.app_context():
		db.create_all()
	app.run()
