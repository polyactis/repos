import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render
from helloworld import model

from pylons.decorators import jsonify
log = logging.getLogger(__name__)
import sys, traceback

class WebhelpController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/template.mako')
		# or, Return a response
		return 'Hello World'

	def fetch(self, id=None):
		helpID = request.params.get('helpID')
		row = model.Stock_250kDB.WebHelp.query.filter_by(short_name=helpID).first()
		if row:
			return row.content
		else:
			return ''
	
	def saveOrUpdate(self, id=None):
		helpID = request.params.get('helpID')
		helpContent = request.body
		row = model.Stock_250kDB.WebHelp.query.filter_by(short_name=helpID).first()
		if row:
			row.content = helpContent
		else:
			row = model.Stock_250kDB.WebHelp(short_name=helpID, content=helpContent)
		try:
			model.db.session.save_or_update(row)
			model.db.session.flush()
			return "1"
		except:
			raise
			return repr(sys.exc_info()) + "\n" + repr(traceback.print_exc())