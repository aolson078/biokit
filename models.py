from Bio import Entrez
from custom_exceptions import InvalidQueryException
from app import create_app, db
from flask_login import UserMixin


# UserMixin gives us premade methods for authentication and user management
class User(UserMixin, db.Model):
	__tablename__ = 'user'
	# for ease of use
	table = __tablename__

	id = db.Column(db.Integer, primary_key=True)
	username = db.Column(db.String(80), unique=True, nullable=False)
	email = db.Column(db.String(120), unique=True, nullable=False)
	password = db.Column(db.String(300), nullable=False, unique=True)

	# array that will hold all currently active (but separate) report for each User. Pickle serializes it for storage
	# report_list = db.Column(db.PickleType, nullable=True)
	# num_documents = 0

	# looks up user in db and adds document to list
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
	nucleotide_id = db.Column(db.String(20), unique=True, nullable=False)
	organism = db.Column(db.String(80), nullable=False)
	gene_info = db.Column(db.String(100), nullable=False)
	allele_info = db.Column(db.String(100))
	sequence_info = db.Column(db.String(100))
	nucleotides = db.Column(db.Text, nullable=False)

	# SAMPLE COMPUTED DATA STORE NOT FINAL PART OF RECORD, (stuff that is calculated w only 1 record)
	# siRNA = db.Column(db.Text, nullable=False)
	# nuc_string_index = db.Column(db.Text, nullable=False)
	# secondary_structure_prediction = db.Column(db.String(100))

	def __repr__(self):
		return '<Organism: %r, NucID: %r, GeneInfo: %r>' % self.organism, self.nucleotide_id, self.gene_info


# the report class represents the final product. It will contain the computed data from the bio processes
class Report(db.Model):
	__tablename__ = 'report'
	id = db.Column(db.Integer, primary_key=True)
	nucleotide_id = db.Column(db.String(20), unique=True, nullable=False)
	organism = db.Column(db.String(80), nullable=False)
	gene_info = db.Column(db.String(100), nullable=False)
	allele_info = db.Column(db.String(100))
	sequence_info = db.Column(db.String(100))
	nucleotides = db.Column(db.Text, nullable=False)

# SAMPLE COMPUTED DATA STORE NOT FINAL PART OF RECORD, (stuff that is calculated w > 1 record)
# phylo_tree = db.Column(db.PickleType nullable=True)


# dot_line_graph = db.Column(db.PickleType nullable=True)


# queries the selected database for term and returns the record with the nucleotide string
def fetch_records(query):
	Entrez.email = "aolson078@gmail.com"
	Entrez.apikey = "4a1d5a80996f3691a67335ded0b79299c708"
	# nucleotide, gene, or protein for db
	database = "nucleotide"
	# number of desired responses per query
	return_max = 3
	# fetch nucleotide record related to term
	IDs = Entrez.read(Entrez.esearch(db=database, term=query, field="Organism", retmax=return_max))["IdList"]
	records = []
	for ID in IDs:
		handle = Entrez.efetch(db=database, id=ID, rettype="fasta", retmod="text")
		record = handle.read()
		records.append(record)
	if len(records) == 0:
		raise InvalidQueryException(query, database)
	return records


# adds and commits new record to database
def add_to_db(record):
	db.session.add(record)
	db.session.commit()


# parses fasta file and creates a database entry for each query result
def parse_records(query):
	records = fetch_records(query)
	for record in records:
		lines = record.split("\n")
		# first line contains the description
		description_line = lines[0]
		# rest of the lines contains the nucleotide sequence
		nucleotides = "".join(lines[1:])

		info = description_line.split(" ")
		nucleotide_id = info[0][1:]  # Remove the '>' character
		organism = " ".join(info[1:3])
		description = " ".join(info[3:]).split(",")

		if len(description) == 3:
			gene_info = description[0]
			allele_info = description[1][1:]
			sequence_info = description[2][1:]
		elif len(description) == 2:
			gene_info = description[0]
			allele_info = "N/A"
			sequence_info = description[1][1:]
		else:
			gene_info = ' '.join(description)
			allele_info = "N/A"
			sequence_info = "N/A"

		new_record = Record(nucleotide_id=nucleotide_id,
		                    organism=organism,
		                    gene_info=gene_info,
		                    allele_info=allele_info,
		                    sequence_info=sequence_info,
		                    nucleotides=nucleotides)

		add_to_db(new_record)


if __name__ == '__main__':
	app = create_app()
	with app.app_context():
		db.create_all()
	app.run()
