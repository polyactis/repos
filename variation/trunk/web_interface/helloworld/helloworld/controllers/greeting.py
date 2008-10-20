import logging

from helloworld.lib.base import *
#from webhelpers.html import literal

log = logging.getLogger(__name__)
import helloworld.model as model

class GreetingController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		c.greeting = request.params.get('greeting', '<b>Welcome</b>')
		#2008-10-05 test db connection
		#row = model.Stock_250kDB.GeneListType.query.first()
		
		c.name = request.params.get('name', 'Visitor')
		return render('/greeting.html')

