import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render
#from helloworld import model

log = logging.getLogger(__name__)

class SnpController(BaseController):

	def index(self, chromosome=None, position=None, call_method_id=None, phenotype_method_id=None):
		"""
		2009-2-19
			snp context
			snp annotation
			snp frequency
			geographic distribution colored by phenotype
			which accession has which allele
			phenotype distribution stratified according to allele
			significant associations (top 1k) in other phenotypes
			
		"""
		# Return a rendered template
		#   return render('/template.mako')
		# or, Return a response
		return 'Hello World'