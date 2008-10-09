import logging

from helloworld.lib.base import *

log = logging.getLogger(__name__)
import helloworld.model as model
class CallinfoController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		c.call_info_ls = []
		return render('/call_info.html')
	
	def call_method(self, id=None):
		query = "select * from view_call"
		if id:
			query += ' where call_method_id=%s'%id
		query += ' order by nativename, stockparent'
		rows = model.db.metadata.bind.execute(query)
		
		c.call_info_ls = []
		i = 0
		for row in rows:
			i += 1
			row.no = i
			c.call_info_ls.append(row)
		return render('/call_info.html')
