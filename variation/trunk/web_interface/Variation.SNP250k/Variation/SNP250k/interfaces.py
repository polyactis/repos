from zope.interface import Interface, Attribute
from zope.app.container.constraints import contains

from zope import schema
from Variation.SNP250k import VariationMessageFactory as _

class PassingData(object):
	pass

class IPhenotypeMethod(Interface):
	"""
	PhenotypeMethod Interface
	"""	
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
from zope.schema.interfaces import IVocabularyFactory
from zope.interface import directlyProvides
from zope.component import getUtility
def MethodIDVocabulary(context):
	locator = getUtility(IDBLocator)
	v = locator.get_phenotype_method_id_ls()
	return SimpleVocabulary.fromItems(v)
directlyProvides(MethodIDVocabulary, IVocabularyFactory)
								
def QCMethodIDVocabulary(context):
	locator = getUtility(IDBLocator)
	v = locator.get_QC_method_id_ls()
	return SimpleVocabulary.fromItems(v)
directlyProvides(QCMethodIDVocabulary, IVocabularyFactory)

def CallMethodIDVocabulary(context):
	locator = getUtility(IDBLocator)
	v = locator.get_call_method_id_ls()
	return SimpleVocabulary.fromItems(v)
directlyProvides(CallMethodIDVocabulary, IVocabularyFactory)

def ResultsMethodTypeIDVocabulary(context):
	locator = getUtility(IDBLocator)
	v = locator.get_results_method_type_id_ls()
	return SimpleVocabulary.fromItems(v)
directlyProvides(ResultsMethodTypeIDVocabulary, IVocabularyFactory)

class IPhenotype(Interface):
	"""
	Phenotype Interface
	"""
	title = schema.TextLine(title=u"Title", 
							required=True)
							
	description = schema.SourceText(title=u"Description", 
								  description=u"A short summary for this phenotype checking out.",
								  required=False)
	#method_id_ls = schema.Choice only picks up simple one-selection widget
	#MultiSelectWidget requires a schema.Choice embedded in a schema.List http://wiki.zope.org/zope3/FAQ2007#how-can-i-use-multiselectwidget-with-formlib
	method_id_ls = schema.List(title=u'Method ID list', required=True,  description=u'Phentoype Method ID in phenotype_method',\
							value_type=schema.Choice(__name__='method_id', title=u'Method ID Choice', vocabulary="Variation.SNP250k.MethodIDVocabulary"))
	
	#method_id_short_name = schema.TextLine(title=u'Method ID Short Name',
	#					 description=u'Phentoype Method ID and Short Name',
	#					 required=True)
	
	phenotype_obj_ls = schema.List(title=u'phenotype object list',
							   description=u"list of phenotype objects, each object has attibute id, short_name, method_description",
							   required=True)
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
	
	text = schema.SourceText(title=u"Text", description=u"More descriptive text.",
						   required=False)

class IDBLocator(Interface):
	"""A utility used to locate appropriate screenings based on search criteria
	"""
	
	def get_phenotype_method_id_ls():
		"""
		return a vocubulary
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

class IQCMethod(Interface):
	"""
	QCMethod Interface
	"""	
	short_name = schema.TextLine(title=u"Short Name",
                            required=True)
	
	method_description = schema.SourceText(title=u"QC Method Description",
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

class IQCOnDirectory(Interface):
	"""
	Interface
	"""
	title = schema.TextLine(title=u"Title",
							required=True)
							
	description = schema.TextLine(title=_(u"Description"),
								  description=_(u"A short summary for this QC."))
	
	QC_method_id = schema.Choice(title=u'QC Method ID',
						 description=u'QC Method Identifier',
						 required=True, vocabulary="Variation.SNP250k.QCMethodIDVocabulary")
	
	input_dir = schema.TextLine(title=u'Input Directory',
						 description=u'Directory Containing Call Files',
						 required=True)
	
	short_name = schema.TextLine(title=u'Short Name',
							   description=u"short name for this QC Method",
							   required=True)
	
	row_id2NA_mismatch_rate = schema.Dict(title=u"QC results. NA_rate and mismatch_rate",
							description=u"dictionary storing results",
							required = True)

class IResults2DB_250k(Interface):
	short_name = schema.TextLine(title=u'Short Name',
							   description=u"short name for this result. Must be unique from previous ones. combining your name, phenotype, data, method is a good one.",
							   required=True)
		
	username = schema.ASCIILine(title=u"Database Username", 
								description=u"This will overwrite/change the setting in DBSetting tab. If your db account doesn't have insert privilege, leave this blank.",
								required=False)

	password = schema.Password(title=u"Database Password", 
								description=u"This will overwrite/change the setting in DBSetting tab.",
								required=False)
	
	phenotype_method_id = schema.Choice(title=u'Phenotype Method ID',
						 description=u'Which Phenotype Used. Pick "no value" if no phenotype used.',
						 required=False, vocabulary="Variation.SNP250k.MethodIDVocabulary")
	"""
	phenotype_method_id = schema.List(title=u'Phenotype Method ID', required=False, \
									description=u'Which Phenotype Used',\
							value_type=schema.Choice(__name__='phenotype_method_id', title=u'Phenotype Method ID Choice', vocabulary="Variation.SNP250k.MethodIDVocabulary"))
	"""
	call_method_id = schema.Choice(title=u'Call Method ID',
						 description=u'From which call method this data is derived from',
						 required=True, vocabulary="Variation.SNP250k.CallMethodIDVocabulary")
	data_description = schema.SourceText(title=u'Data Description',
							   description=u"Describe how your data is derived from that call method. like non-redundant set, 1st 96, etc.",
							   required=True)
	results_method_type_id = schema.Choice(title=u'Results Method Type ID', required=False, \
									description=u'Which type of results. Pick "no value" if no appropriate type. Fill in the blank right below.',\
									vocabulary="Variation.SNP250k.ResultsMethodTypeIDVocabulary")
	results_method_type_short_name = schema.TextLine(title=u'New Results Method Type Short Name', required=False, \
									description=u'To create a new results method type if the above does not have your results type.')
	
	method_description = schema.SourceText(title=u'Method Description',
							   description=u"Describe your method and what type of score, association (-log or not), recombination etc.",
							   required=True)
	comment = schema.SourceText(title=u'Further Comment',
							   description=u"Anything worth other people to know?",
							   required=False)
	input_fname = schema.Bytes(title=u"Input File", 
								  description=u"File containing Genome-Wide Results. tab or comma delimited 3-column (chromosome, position, score) or 4-column (chromosome, start_position, end_position, score).",
								  required=True)
	commit_type = schema.Choice(title=u'Commit Type',
						description=u'Which type of database commit action you want.',
						required=True, vocabulary=SimpleVocabulary.fromItems([("No commit. Test Run!", 0),\
																			("Commit this transaction.", 1)]))
						#required=True, vocabulary=SimpleVocabulary.fromItems([("No commit this transaction. Leave it to the next transaction. Save Time!", 0),\
						#													("Commit this and all previous transactions.", 1),\
						#													("Commit all previous (excluding the current page) transactions. Fill in the required field in the current page with whatever. They won't go into database.", 2),\
						#													("Rollback all previous transactions. This = Undo!", 3)]))