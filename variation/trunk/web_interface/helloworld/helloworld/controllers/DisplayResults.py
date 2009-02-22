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
import gviz_api	#2009-2-1 google visualization python api to export data into json-related formats
import numpy, datetime

class DisplayresultsController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/template.mako')
		# or, Return a response
		return 'Hello World'
	
	@staticmethod
	def getCallMethodLs():
		"""
		2009-1-30
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
		2009-1-30
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
		2009-1-30
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
		2009-1-30
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
		2009-1-30
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
	
	def getCallMethodEigenValue(self):
		"""
		2009-2-2
			get an eigen value given a call_method_id
		"""
		call_method_id = request.params.getone('call_method_id')
		which_eigen_value = request.params.getone('which_eigen_value')
		rows = model.Stock_250kDB.CallMethodEigenValues.query.filter_by(call_method_id=call_method_id).filter_by(which_eigen=which_eigen).all()
	
	def getCallInfoData(self, call_method_id, phenotype_method_id):
		"""
		2009-2-2
			get ecotype id, name, latitude, longitude, pc1, pc2, phenotype for all call_infos in one call_method_id
		"""
		#1st get call_info_id 2 PC12 mapping
		rows = model.Stock_250kDB.CallInfo.query.filter_by(method_id=call_method_id).all()
		call_info_id2pc12 = {}
		for row in rows:
			call_info_id2pc12[row.id] = row.pc_value_ls[:2]
		
		
		#2nd check if model.pheno_data exists
		pheno_data = getattr(model, 'pheno_data', None)
		if pheno_data is None:
			from variation.src.OutputPhenotype import OutputPhenotype
			pheno_data = OutputPhenotype.getPhenotypeData(model.db.metadata.bind, phenotype_avg_table=model.Stock_250kDB.PhenotypeAvg.table.name,\
														phenotype_method_table=model.Stock_250kDB.PhenotypeMethod.table.name)
			model.pheno_data = pheno_data
		#3rd finally construct the full data and turn it into json
		column_name_type_ls = [("date", ("date", "Date")), ("ecotypeid", ("number", "Ecotype ID")), ("label", ("string","ID Name Phenotype")), \
							("lat",("number", "Latitude")), ("lon",("number", "Longitude")), ("name", ("string", "Native Name")), \
							("pc1", ("number","PC1")), ("pc2", ("number", "PC2")), ("phenotype", ("number", "Phenotype"))]
		
		description = dict(column_name_type_ls)
		rows = model.db.metadata.bind.execute("select * from view_call where call_method_id=%s"%(call_method_id))
		return_ls = []
		for row in rows:
			pc12 = call_info_id2pc12.get(row.call_info_id)
			if pc12:
				call_info_pc1, call_info_pc2 = pc12[:2]
				pc1 = call_info_pc1.pc_value
				pc2 = call_info_pc2.pc_value
			else:
				pc1 = pc2 = 0
			pheno_row_index = pheno_data.row_id2row_index.get(row.ecotype_id)
			pheno_col_index = pheno_data.col_id2col_index.get(int(phenotype_method_id))
			if pheno_row_index is not None and pheno_col_index is not None:
				phenotype_value = pheno_data.data_matrix[pheno_row_index][pheno_col_index]
				if phenotype_value == 'NA':
					phenotype_value = -0.01
			else:
				phenotype_value = -0.01	#numpy.nan can't be recognized by ToJSon()
			label = 'ID:%s. Name=%s. Phenotype=%s.'%(row.ecotype_id, row.nativename, phenotype_value)
			return_ls.append(dict(ecotypeid=row.ecotype_id, name=row.nativename, lat=row.latitude, lon=row.longitude,\
								pc1=pc1, pc2=pc2, phenotype=phenotype_value, label=label, date=datetime.date(2009,2,3)))
		
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		#column_ls.sort()	#['date', 'ecotypeid', 'label', 'lat', 'lon', 'name', 'pc1', 'pc2', 'phenotype']
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		return json_result
	
	def fetchCallInfoData(self):
		"""
		2009-2-2
			get a complete catalog of call information given a call_method_id
			including, ecotypeid, name, latitude, longitude, phenotype value, PC1, PC2
		"""
		call_method_id = request.params.getone('call_method_id')
		phenotype_method_id = request.params.getone('phenotype_method_id')
		
		key_tup = (call_method_id, phenotype_method_id)
		CallPhenotypeMethod2call_info_jsdata = getattr(model, 'CallPhenotypeMethod2call_info_data', None)
		if CallPhenotypeMethod2call_info_jsdata is None:
			model.CallPhenotypeMethod2call_info_jsdata = {}
			CallPhenotypeMethod2call_info_jsdata = model.CallPhenotypeMethod2call_info_jsdata
		
		call_info_jsdata = CallPhenotypeMethod2call_info_jsdata.get(key_tup)
		if call_info_jsdata is None:
			call_info_jsdata = self.getCallInfoData(call_method_id, phenotype_method_id)
		return call_info_jsdata
	
	def getPhenotypeHistImage(self, id=None, img_type='hist_thumb'):
		phenotype_method_id = request.params.get('phenotype_method_id', id)
		img_type = request.params.get('img_type', img_type)
		
		if phenotype_method_id is None:
			response.headers['Content-type'] = 'text/plain'
			return 'No image'
		phenotype_method_id = int(phenotype_method_id)
		
		img_type = request.params.get('img_type', img_type)
		if img_type.find('svg')!=-1:
			response.headers['Content-type'] = 'image/svg'
		else:
			response.headers['Content-type'] = 'image/png'
		
		get_img_data_success = 0
		if phenotype_method_id:
			row = model.Stock_250kDB.PhenotypeHistPlots.query.filter_by(phenotype_method_id=phenotype_method_id).first()
			if row:
				img_data = getattr(row, img_type, None)
				if img_data:
					img_data = img_data.__str__()	#2008-12-26	img_data is a buffer. weird!
					get_img_data_success = 1
				else:
					get_img_data_success = 0
		if not get_img_data_success:
			response.headers['Content-type'] = 'text/plain'
			img_data = 'No Image'
		return img_data
	
	def getCallPhenotypeQQImage(self, id=None, img_type='qq_thumb'):
		phenotype_method_id = request.params.get('phenotype_method_id', id)
		call_method_id = request.params.get('call_method_id', 29)
		img_type = request.params.get('img_type', img_type)
		
		if phenotype_method_id is None:
			response.headers['Content-type'] = 'text/plain'
			return 'No image'
		phenotype_method_id = int(phenotype_method_id)
		call_method_id = int(call_method_id)
		
		img_type = request.params.get('img_type', img_type)
		if img_type.find('svg')!=-1:
			response.headers['Content-type'] = 'image/svg'
		else:
			response.headers['Content-type'] = 'image/png'
		
		get_img_data_success = 0
		if call_method_id and phenotype_method_id:
			row = model.Stock_250kDB.CallPhenotypeQQPlots.query.filter_by(phenotype_method_id=phenotype_method_id).\
						filter_by(phenotype_method_id=phenotype_method_id).first()
			if row:
				img_data = getattr(row, img_type, None)
				if img_data:
					img_data = img_data.__str__()	#2008-12-26	img_data is a buffer. weird!
					get_img_data_success = 1
				else:
					get_img_data_success = 0
		if not get_img_data_success:
			response.headers['Content-type'] = 'text/plain'
			img_data = 'No Image'
		return img_data