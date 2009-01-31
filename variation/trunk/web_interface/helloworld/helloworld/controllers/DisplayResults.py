import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render
#from helloworld import model

log = logging.getLogger(__name__)

from pymodule import PassingData
from formencode import htmlfill
from pylons.decorators import jsonify	#2008-12-30
from helloworld.lib.base import *

class DisplayresultsController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/template.mako')
		# or, Return a response
		return 'Hello World'
	
	@staticmethod
	def getCallMethodLs():
		"""
		2008-12-30
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name
		list_info  = h.getCallMethodInfo(affiliated_table_name=affiliated_table_name)
		call_method_ls = []
		for i in range(len(list_info.id_ls)):
			id = list_info.id_ls[i]
			label = list_info.label_ls[i]
			call_method_ls.append([id, label])
		return call_method_ls
	
	@staticmethod
	def getPhenotypeMethodLs(call_method_id):
		"""
		2008-12-30
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name
		extra_condition = 's.call_method_id=%s'%(call_method_id)
		phenotype_info  = h.getPhenotypeInfo(affiliated_table_name=affiliated_table_name, extra_condition=extra_condition)
		phenotype_method_ls = []
		for i in range(len(phenotype_info.phenotype_method_id_ls)):
			phenotype_method_id = phenotype_info.phenotype_method_id_ls[i]
			phenotype_method_label = phenotype_info.phenotype_method_label_ls[i]
			phenotype_method_ls.append([phenotype_method_id, phenotype_method_label])
		return phenotype_method_ls
	
	@jsonify
	def getPhenotypeMethodLsJson(self):
		result = {
				'options': [
						dict(id=value, value=id) for id, value in self.getPhenotypeMethodLs(request.params.getone('call_method_id'))
						]
				}
		#result['options'].append({'id': u'[At the end]', 'value': u''})
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return result
	
	@staticmethod
	def getAnalysisMethodLs(call_method_id, phenotype_method_id):
		"""
		2008-12-30
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name	#alias is 's'
		extra_condition = 's.call_method_id=%s and s.phenotype_method_id=%s'%\
			(call_method_id, phenotype_method_id)
		list_info = h.getAnalysisMethodInfo(affiliated_table_name, extra_condition=extra_condition)
		
		analysis_method_ls = []
		
		for i in range(len(list_info.id_ls)):
			id = list_info.id_ls[i]
			label = list_info.label_ls[i]
			analysis_method_ls.append([id, label])
		return analysis_method_ls
	
	@jsonify
	def getAnalysisMethodLsJson(self):
		call_method_id = request.params.getone('call_method_id')
		phenotype_method_id = request.params.getone('phenotype_method_id')
		result = {
				'options': [
						dict(id=value, value=id) for id, value in self.getAnalysisMethodLs(call_method_id, phenotype_method_id)
						]
				}
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return result
	
	def form(self, id=None):
		"""
		2008-12-30
		"""
		if id is None:
			id = request.params.get('id', 1)
		defaults = {'call_method_id': id,
				"remove_old_plots":"on"}
		c.call_method_ls = self.getCallMethodLs()
		c.call_method_ls.insert(0, [0, u'Please Choose ...'])
		html = render('/display_results_form.html')
		return htmlfill.render(html, defaults)
	
	@jsonify
	def fetchOne(self, id=None):
		"""
		"""
		if id is None:
			id = request.params.get('id', None)
		if id is None:
			id = request.params.get('results_id', None)
		
		call_method_id = request.params.getone('call_method_id')
		phenotype_method_id = request.params.getone('phenotype_method_id')
		analysis_method_id = request.params.getone('analysis_method_id')
		rm = model.Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).\
			filter_by(phenotype_method_id=phenotype_method_id).filter_by(analysis_method_id=analysis_method_id).first()
		result_id = rm.id
		
		param_data = PassingData(min_MAF=0.1)
		
		genome_wide_result = h.GeneListRankTest.getResultMethodContent(rm, min_MAF=0.1, pdata=param_data)
		
		return_ls = []
		import gviz_api
		description = {"position": ("number", "Position"),"value": ("number", "-log Pvalue")}
		chr2data = {}
		for i in range(10000):
			data_obj = genome_wide_result.get_data_obj_at_given_rank(i+1)
			if data_obj.chromosome not in chr2data:
				chr2data[data_obj.chromosome] = []
			chr2data[data_obj.chromosome].append(dict(position=data_obj.position, value=data_obj.value))
		
		json_result = {}
		for chr in chr2data:
			data_table = gviz_api.DataTable(description)
			data_table.LoadData(chr2data[chr])
			json_result[chr] = data_table.ToJSon(columns_order=("position", "value"))
												#,order_by="position")
		"""
		for i in range(10000):
			data_obj = genome_wide_result.get_data_obj_at_given_rank(i+1)
			return_ls.append(dict(chromosome=data_obj.chromosome, position=data_obj.position, value=data_obj.value))
		
		result = {
				'gwr': return_ls
				}
		"""
		return json_result
		