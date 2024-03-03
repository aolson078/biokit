from Bio import Entrez
from custom_exceptions import InvalidQueryException
from app import db
from flask_login import UserMixin


# UserMixin gives us premade methods for authentication and user management
class User(UserMixin, db.Model):
	__tablename__ = 'user'

	id = db.Column(db.Integer, primary_key=True)
	username = db.Column(db.String(80), unique=True, nullable=False)
	email = db.Column(db.String(120), unique=True, nullable=False)
	password = db.Column(db.String(300), nullable=False, unique=True)

	def __repr__(self):
		return '<User: %r>' % self.username


# Queries the selected database for term and returns the record with the nucleotide string
def fetch_record(query):
	Entrez.email = "aolson078@gmail.com"
	Entrez.apikey = "4a1d5a80996f3691a67335ded0b79299c708"
	# nucleotide, gene, or protein for db
	database = "nucleotide"
	# Fetch nucleotide record related to term
	IDs = Entrez.read(Entrez.esearch(db=database, term=query, field="Organism", retmax=3))["IdList"]
	records = []
	for ID in IDs:
		handle = Entrez.efetch(db=database, id=ID, rettype="fasta", retmod="text")
		record = handle.read()
		records.append(record)
	if len(records) == 0:
		raise InvalidQueryException(query, database)
	return records
