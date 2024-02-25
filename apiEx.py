from Bio import Entrez
from custom_exceptions import InvalidQueryException
from functions import parse_fasta

Entrez.email = "aolson078@gmail.com"
Entrez.apikey = "4a1d5a80996f3691a67335ded0b79299c708"

database = {1: "nucleotide", 2: "gene", 3: "protein"}


def get_database_choice():
	# Prompt the user for database choice
	print("Choose a database:")
	print("1. Nucleotide")
	print("2. Gene")
	print("3. Protein")
	db_choice = input("Enter the number of the database you wish to query: ")
	return database.get(int(db_choice))


# Queries the selected database for term and returns the record with the nucleotide string
def fetch_record():
	term = input("Please enter term to query (Ex: Homo sapien): ")
	db = get_database_choice()
	# Fetch nucleotide record related to term
	IDs = Entrez.read(Entrez.esearch(db=db, term=term, field="Organism", retmax=1))["IdList"]
	records = []
	for ID in IDs:
		handle = Entrez.efetch(db=db, id=ID, rettype="fasta", retmod="text")
		record = handle.read()
		records.append(record)
	if len(records) == 0:
		raise InvalidQueryException(term, db)
	return records


results = parse_fasta(fetch_record())
print(results[0][0])
print()
print(results[1][0])
