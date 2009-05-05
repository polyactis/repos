import logging

from helloworld.lib.base import *

log = logging.getLogger(__name__)
from HelpOtherControllers import HelpothercontrollersController as hc

class MafvsscoreplotController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		MAFVsScorePlot = model.Stock_250kDB.MAFVsScorePlot
		ResultsMethod = model.Stock_250kDB.ResultsMethod
		rows = model.db.metadata.bind.execute("select r.call_method_id, count(m.id) as count from %s r, %s m where m.results_method_id=r.id\
			group by r.call_method_id"%(ResultsMethod.table.name, MAFVsScorePlot.table.name))
		c.maf_vs_score_summary_ls = []
		for row in rows:
			call_method = model.Stock_250kDB.CallMethod.get(row.call_method_id)
			maf_vs_score_summary = PassingData(call_method_id=row.call_method_id, call_method_short_name=call_method.short_name, no_of_plots=row.count)
			c.maf_vs_score_summary_ls.append(maf_vs_score_summary)
		return render('/maf_vs_score.html')
	
	def call_method(self, id=None):
		if id is None:
			id = request.params.get('id', 17)
		
		MAFVsScorePlot = model.Stock_250kDB.MAFVsScorePlot
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name
		extra_condition = 's.call_method_id=%s and s.id=m.results_method_id'%id
		extra_tables = '%s m'%model.Stock_250kDB.MAFVsScorePlot.table.name
		c.phenotype_info  = hc.getPhenotypeInfo(affiliated_table_name, extra_condition, extra_tables)
		c.analysis_info = hc.getAnalysisMethodInfo(affiliated_table_name, extra_condition, extra_tables)
		
		rows = MAFVsScorePlot.query.filter(MAFVsScorePlot.results_method.has(call_method_id=id))
		data_matrix = numpy.zeros([len(c.phenotype_info.phenotype_method_id_ls), len(c.analysis_info.id_ls)], numpy.int)
		data_matrix[:] = -1
		counter = 0
		for row in rows:
			row_index = c.phenotype_info.phenotype_method_id2index[row.results_method.phenotype_method_id]
			col_index = c.analysis_info.id2index[row.results_method.analysis_method_id]
			data_matrix[row_index, col_index] = row.id
			counter +=1
		c.counter = counter
		c.data_matrix = data_matrix
		c.call_method_id = id
		return render('/maf_vs_score_call_method.html')
	
	def getImage(self, id=None, img_type='png_data'):
		img_type = request.params.get('img_type', img_type)
		if img_type.find('svg')!=-1:
			response.headers['Content-type'] = 'image/svg'
		else:
			response.headers['Content-type'] = 'image/png'
		if id is None:
			id = request.params.get('id', None)
		get_img_data_success = 0
		if id:
			entry = model.Stock_250kDB.MAFVsScorePlot.get(id)
			if entry:
				img_data = getattr(entry, img_type, None)
				if img_data:
					get_img_data_success = 1
					img_data = img_data.__str__()	#2008-12-26	img_data becomes a buffer now. weird!
				else:
					get_img_data_success = 0
		if not get_img_data_success:
			response.headers['Content-type'] = 'text/plain'
			img_data = 'No Image'
		return img_data