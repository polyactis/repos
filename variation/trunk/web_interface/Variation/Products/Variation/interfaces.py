from zope.interface import Interface, Attribute
from zope.app.container.constraints import contains

from zope import schema

class IPhenotypeMethod(Interface):
	"""
	PhenotypeMethod Interface
	"""
	id = schema.Int(title=u"Phenotype Method identifier",
							  description=u"A unique id for this phenotype_method",
							  required=True,
							  readonly=True)
	
	short_name = schema.TextLine(title=u"Short Name",
                            required=True)
	
	method_description = schema.SourceText(title=u"Phenotype Method Description",
						   readonly=True)
	data_description = schema.SourceText(title=u"Data Description",
						   readonly=True)
	comment = schema.SourceText(title=u"Comment",
						   readonly=True)
	created_by = schema.TextLine(title=u"Created By",
                            readonly=True)
	updated_by = schema.TextLine(title=u"Updated By",
                            readonly=True)
	date_created = schema.Datetime(title=u"Date Created",
							readonly=True)
	date_updated = schema.Datetime(title=u"Date Updated",
							readonly=True)

class IPhenotypeAvg(Interface):
	"""
	PhenotypeAvg Interface
	"""
	ecotype_id = schema.Int(title=u"Ecotype id",
							  description=u"A unique id for the ecotype",
							  required=True,
							  readonly=True)
	value = schema.Float(title=u"Phenotype Value",
							  description=u"",
							  readonly=True)
	stdev = schema.Float(title=u"Phenotype Standard Deviation",
							  description=u"",
							  readonly=True)
	sample_size = schema.Int(title=u"Phenotype Standard Deviation",
							  description=u"",
							  readonly=True)
	method_id = schema.Int(title=u"Phenotype Method ID",
							  description=u"",
							  required=True,
							  readonly=True)

from zope.schema.vocabulary import SimpleVocabulary
from zope.component import getUtility
def MethodIDVocabulary(context):
	locator = getUtility(IPhenotypeLocator)
	v = locator.get_phenotype_method_id_ls()
	return SimpleVocabulary.fromItems(v)


class IPhenotype(Interface):
	"""
	Phenotype Interface
	"""
	title = schema.TextLine(title=u"Title", 
							required=True)
							
	description = schema.SourceText(title=u"Description", 
								  description=u"A short summary for this phenotype checking out.")
	
	method_id_ls = schema.Choice(title=u'Method ID',
						 description=u'Phentoype Method ID in phenotype_method',
						 required=True, values=['8499','8257'])
	
	#method_id_short_name = schema.TextLine(title=u'Method ID Short Name',
	#					 description=u'Phentoype Method ID and Short Name',
	#					 required=True)
	
	short_name_ls = schema.List(title=u'Short Name',
							   description=u"short name for this phenotype",
							   required=True)
	
	method_description_ls = schema.List(title=u"Method Description",
                             description=u"Descriptive text about this phenotype.")
	
	data_matrix = schema.List(title=u"Data Matrix",
                             description=u"Data Matrix. each row is ecotype_id, nativename, value.")
	"""
	data_description = schema.SourceText(title=u"Data Description",
						   readonly=True)
	comment = schema.SourceText(title=u"Comment",
						   readonly=True)
	created_by = schema.TextLine(title=u"Created By",
                            readonly=True)
	updated_by = schema.TextLine(title=u"Updated By",
                            readonly=True)
	date_created = schema.Datetime(title=u"Date Created",
							readonly=True)
	date_updated = schema.Datetime(title=u"Date Updated",
							readonly=True)
	"""
	
class IResults(Interface):
	"""
	Results Interface
	"""

class ICelData(Interface):
	"""
	CelData Interface
	"""
class IGenotype(Interface):
	"""
	Genotype Interface
	"""



class IVariationFolder(Interface):
	"""A folder containing phenotype, results, cel data, genotype calls ...
	"""
	contains('Variation.interfaces.IPhenotype', 
			 'Variation.interfaces.IResults', 
			 'Variation.interfaces.ICelData', 
			 'Variation.interfaces.IGenotype',)
	
	title = schema.TextLine(title=u"Title", 
							required=True)
							
	description = schema.TextLine(title=u"Description", 
								  description=u"A short summary of this folder")

class IPhenotypeLocator(Interface):
	"""A utility used to locate appropriate screenings based on search criteria
	"""
	
	def films_at_cinema(cinema, from_date, to_date):
		"""Return a list of all films screening at the particular ICinema
		between the specified dates.
		
		Returns a list of dictionaries with keys 'film_code', 'url', 'title' 
		and 'summary'.
		"""
		
	def cinemas_for_film(film, from_date, to_date):
		"""Return a list of all cinemas screening the given film between the
		specified dates.
		
		Returns a list of dictionaries with keys 'cinema_code', 'url', 'name' 
		and 'address'.
		"""
		
	def screenings(film, cinema, from_date, to_date):
		"""Return all screenings of the given film, at the given cinema,
		between the given dates
		
		Returns a list of IScreening objects.
		"""
		
	def screening_by_id(screening_id):
		"""Get an IScreening from a screening id
		"""

class PhenotypeError(Exception):
	"""Exception raised if there is an error making a reservation
	"""

	def __init__(self, message):
		Exception.__init__(self, message)
		self.error_message = message

class IStockDatabaseSettings(Interface):
	"""Database connection settings.
	"""
	
	drivername = schema.ASCIILine(title=u"Driver name", 
								  description=u"The database driver name", 
								  default='mysql', 
								  required=True)

	hostname = schema.ASCIILine(title=u"Host name", 
								description=u"The database host name", 
								default='localhost', 
								required=True)
								
	port = schema.Int(title=u"Port number", 
					  description=u"The database port number. Leave blank to use the default.", 
					  required=False)
								
	username = schema.ASCIILine(title=u"User name", 
								description=u"The database user name", 
								required=True)

	password = schema.Password(title=u"Password", 
								description=u"The database password", 
								required=False)
								
	database = schema.ASCIILine(title=u"Database name", 
								description=u"The name of the database on this server", 
								required=True)