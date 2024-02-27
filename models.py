from Bio import Entrez
from custom_exceptions import InvalidQueryException


# Queries the selected database for term and returns the record with the nucleotide string
def fetch_record(query):
	Entrez.email = "aolson078@gmail.com"
	Entrez.apikey = "4a1d5a80996f3691a67335ded0b79299c708"
	# nucleotide, gene, or protein for db
	db = "nucleotide"
	# Fetch nucleotide record related to term
	IDs = Entrez.read(Entrez.esearch(db=db, term=query, field="Organism", retmax=3))["IdList"]
	records = []
	for ID in IDs:
		handle = Entrez.efetch(db=db, id=ID, rettype="fasta", retmod="text")
		record = handle.read()
		records.append(record)
	if len(records) == 0:
		raise InvalidQueryException(query, db)
	return records
