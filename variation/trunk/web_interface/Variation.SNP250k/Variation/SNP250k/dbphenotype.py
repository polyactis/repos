from zope.interface import implements
from zope.component import getUtility

from Variation.SNP250k.interfaces import IPhenotypeMethod, IPhenotypeAvg,  IQCMethod
from Variation.SNP250k.interfaces import IPhenotypeLocator
#from optilux.cinemacontent.interfaces import ITicketReservations
#from optilux.cinemacontent.interfaces import ReservationError

#from optilux.cinemacontent.screening import Screening
#from Products.Variation import VariationMessageFactory as _

import sqlalchemy as sql
from collective.lead.interfaces import IDatabase

class PhenotypeMethod(object):
	"""PhenotypeMethod for ORM mapping
	"""
	
	implements(IPhenotypeMethod)
	
	id = None
	short_name = None
	method_description = u""
	data_description = u""
	comment = u""
	created_by = None
	updated_by = None
	date_created = None
	date_updated = None
	
	def __init__(self, id, short_name, method_description, data_description, comment, created_by, updated_by, date_created, date_updated):
		self.id = id
		self.short_name = short_name
		self.method_description = method_description
		self.data_description = data_description
		self.comment = comment
		self.created_by = created_by
		self.updated_by = updated_by
		self.date_created = date_created
		self.date_updated = date_updated

class PhenotypeAvg(object):
	"""PhenotypeAvg for ORM mapping
	"""
	
	implements(IPhenotypeAvg)
	id = None
	ecotype_id = None
	value = None
	stdev = None
	sample_size = None
	method_id = None
	
	def __init__(self, id, ecotype_id, value, stdev, sample_size, method_id):
		self.id = id
		self.ecotype_id = ecotype_id
		self.value = value
		self.stdev = stdev
		self.sample_size = sample_size
		self.method_id = method_id

class QCMethod(object):
	"""QCMethod for ORM mapping
	"""
	
	implements(IQCMethod)
	
	short_name = None
	method_description = u""
	data_description = u""
	comment = u""
	created_by = None
	updated_by = None
	date_created = None
	date_updated = None
	
	def __init__(self, id, short_name, method_description, data_description, comment, created_by, updated_by, date_created, date_updated):
		self.id = id
		self.short_name = short_name
		self.method_description = method_description
		self.data_description = data_description
		self.comment = comment
		self.created_by = created_by
		self.updated_by = updated_by
		self.date_created = date_created
		self.date_updated = date_updated

#class TicketReservations(object):
#	"""Make reservations in the reservations database
#	"""
#	
#	implements(ITicketReservations)
#	
#	def __call__(self, reservation):
#		"""Make a reservation
#		"""
#		
#		db = getUtility(IDatabase, name='optilux.reservations')
#		session = db.session
#		
#		# Make sure there are still seats available
#		screening = reservation.screening
#		session.refresh(screening)
#		
#		if screening.remaining_tickets <= 0:
#			raise ReservationError(_(u"This screening is sold out!"))
#		elif screening.remaining_tickets < reservation.num_tickets:
#			raise ReservationError(_(u"Not enough tickets remaining!"))
#		
#		# Otherwise, we can save the reservation
#		screening.remaining_tickets -= reservation.num_tickets
#		session.update(screening)
#		session.save(reservation)
#		session.flush()

class PhenotypeLocator(object):
	"""Find phenotype in db
	"""
	
	implements(IPhenotypeLocator)
	
	def get_phenotype_method_id_ls(self):
		"""Return a list of all films showing at the particular ICinema
		between the specified dates.
		
		Returns a list of dictionaries with keys 'film_code', 'url', 'title' 
		and 'summary'.
		"""
		
		# Set up and issue the query, making sure we use the same transaction
		# context as SQLAlchemy's session. 
		
		db = getUtility(IDatabase, name='variation.stockdatabase')
		connection = db.connection
				
		statement = sql.select([PhenotypeMethod.c.id, PhenotypeMethod.c.short_name],
							   distinct=True)
		
		results = connection.execute(statement).fetchall()
		
		# Now use the catalog to find films for the returned film codes
		#film_codes = [row['film_code'] for row in results]
		
		# We don't have a context to acquire from here, but we can get
		# it from zope.app.component.hooks.getSite(), which returns the
		# 'site' (in the Component Architecture sense) set during traversal.
		# Since we know that the Plone site root is a "local site" we can
		# be sure that either this, or a sub-site, will be returned when
		# this function is invoked from somewhere within the Plone site.

		#site = getSite()
		#catalog = getToolByName(site, 'portal_catalog')
		"""
		return [ dict(film_code=film.film_code,
					  url=film.getURL(),
					  title=film.Title,
					  summary=film.Description,)
				 for film in 
					catalog(object_provides=IFilm.__identifier__,
							film_code=film_codes,
							sort_on='sortable_title')
			   ]
		"""
		vocabulary = [('%s %s'%(row[0], row[1]), row[0]) for row in results]	#(token, value)
		return vocabulary

	def get_QC_method_id_ls(self):
		"""
		2008-04-23
		"""
		
		# Set up and issue the query, making sure we use the same transaction
		# context as SQLAlchemy's session. 
		
		db = getUtility(IDatabase, name='variation.stockdatabase')
		connection = db.connection
				
		statement = sql.select([QCMethod.c.id, QCMethod.c.short_name],
							   distinct=True, order_by=[QCMethod.c.id])
		
		results = connection.execute(statement).fetchall()
		vocabulary = [('%s %s'%(row[0], row[1]), row[0]) for row in results]	#(token, value)
		return vocabulary
	
	def get_short_name_description_ls(self, method_id_ls):
		"""
		2008-04-18
			method_id_ls is only one long integer (not a list, zope.schema.Choice)
		2008-03-28
		"""
		short_name_ls = []
		method_description_ls = []
		db = getUtility(IDatabase, name='variation.stockdatabase')
		connection = db.connection
		session = db.session
		if type(method_id_ls)!=list:
			method_id_ls = [method_id_ls]
		method_id_ls.sort()
		for method_id in method_id_ls:
			#phenotype_method = session.query(PhenotypeMethod).get(method_id)
			statement = sql.select([PhenotypeMethod.c.short_name, PhenotypeMethod.c.method_description],
							   sql.and_(
									PhenotypeMethod.c.id == method_id),
							   distinct=True)
			results = connection.execute(statement).fetchall()
			#cinema_codes = [row['cinema_code'] for row in results]
			short_name_ls += [row['short_name'] for row in results]
			method_description_ls += [row['method_description'] for row in results]
		
		return short_name_ls, method_description_ls
	
	"""
	def cinemas_for_film(self, film, from_date, to_date):
		db = getUtility(IDatabase, name='optilux.reservations')
		connection = db.connection
		
		statement = sql.select([Screening.c.cinema_code],
							   sql.and_(
									Screening.c.film_code == film.film_code,
									Screening.c.show_time.between(from_date, to_date)
							   ),
							   distinct=True)
		
		results = connection.execute(statement).fetchall()
		
		cinema_codes = [row['cinema_code'] for row in results]
		
		site = getSite()
		catalog = getToolByName(site, 'portal_catalog')
		
		return [ dict(cinema_code=cinema.cinema_code,
					  url=cinema.getURL(),
					  name=cinema.Title,
					  address=cinema.Description,)
				 for cinema in 
					catalog(object_provides=ICinema.__identifier__,
							cinema_code=cinema_codes,
							sort_on='sortable_title')
			   ]
		
	def screenings(self, film, cinema, from_date, to_date):
		db = getUtility(IDatabase, name='optilux.reservations')
		session = db.session
		
		screenings = session.query(Screening).select(sql.and_(
														Screening.c.film_code==film.film_code,
														Screening.c.cinema_code==cinema.cinema_code,
														Screening.c.show_time.between(from_date, to_date)
													 ), order_by=[Screening.c.show_time])

		# Now set the 'film' and 'cinema' properties
		
		for screening in screenings:
			screening.film = film
			screening.cinema = cinema
		
		return screenings
		
	def screening_by_id(self, screening_id):
		db = getUtility(IDatabase, name='optilux.reservations')
		session = db.session
		
		screening = session.query(Screening).get(screening_id)
		
		site = getSite()
		catalog = getToolByName(site, 'portal_catalog')
		
		# Set cinema and film. We assume these are in the portal and that 
		# there is a one-to-one mapping between cinema/film codes and 
		# cinemas/films.

		screening.cinema = catalog(cinema_code=screening.cinema_code)[0].getObject()
		screening.film = catalog(film_code=screening.film_code)[0].getObject()
		
		return screening
	"""
