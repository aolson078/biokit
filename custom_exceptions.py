class InvalidQueryException(Exception):
	"""Exception raised for errors querying the database"""

	def __init__(self, term, db):
		self.term = term
		self.db = db
		self.message = f"No results found for \"{term}\" in {db} database"
		super().__init__(f"Invalid query: {term} in {db}")


class NoRecordsError(Exception):
	"""Exception raised when user tries to create a report without any records in the db"""

	def __init__(self, term):
		self.message = "No records in database"
		super().__init__("No records in database")
