class InvalidQueryException(Exception):
	"""Exception raised for errors querying the database"""

	def __init__(self, term, db):
		self.term = term
		self.db = db
		self.message = f"No results found for \"{term}\" in {db} database"
		super().__init__(f"Invalid query: {term} in {db}")
