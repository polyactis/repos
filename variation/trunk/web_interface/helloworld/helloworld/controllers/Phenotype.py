import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render, config, h
from helloworld import model
from pylons.decorators import jsonify
import simplejson, gviz_api, datetime, numpy
log = logging.getLogger(__name__)
from HelpOtherControllers import HelpothercontrollersController as hc

class PhenotypeController(BaseController):

	def index(self, id=None):
		"""
		2009-4-6
			display a summary page given a phenotype id
		"""
		# Return a rendered template
		#   return render('/template.mako')
		# or, Return a response
		c.phenotype_method_id = request.params.get('phenotype_method_id', id)
		if c.phenotype_method_id:
			c.phenotype_method_id = int(c.phenotype_method_id)
		pm = model.Stock_250kDB.PhenotypeMethod.get(c.phenotype_method_id)
		c.phenotype_method_short_name = pm.short_name
		c.phenotype_method_description = pm.method_description
		
		c.call_method_id = request.params.get('call_method_id', int(config['app_conf']['published_call_method_id']))
		c.callInfoURL = h.url_for(controller='DisplayResults', action='fetchCallInfoData', id=None,\
							phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		c.phenotypeHistImageURL = h.url_for(controller='DisplayResults', action='getPhenotypeHistImage', id=None, \
										phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		c.callPhenotypeQQImageURL = h.url_for(controller='DisplayResults', action='getCallPhenotypeQQImage', id=None,\
											phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		c.phenotypeHistogramDataURL = h.url_for(controller='Phenotype', action='getPhenotypeHistogramData', id=c.phenotype_method_id)
		return render('/OnePhenotype.html')
	
	@jsonify
	def getPhenotypeMethodLs(self, call_method_id=None, affiliated_table_name=None):
		"""
		2009-4-6
			copied from DisplayResults.py
		"""
		#call_method_id = request.params.getone('call_method_id')
		#affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name
		#extra_condition = 's.call_method_id=%s'%(call_method_id)
		
		extra_condition = ''
		phenotype_info  = hc.getPhenotypeInfo(affiliated_table_name=affiliated_table_name, extra_condition=extra_condition)
		phenotype_method_ls = []
		for i in range(len(phenotype_info.phenotype_method_id_ls)):
			phenotype_method_id = str(phenotype_info.phenotype_method_id_ls[i])
			phenotype_method_label = phenotype_info.phenotype_method_label_ls[i]
			phenotype_method_ls.append([phenotype_method_id, phenotype_method_label])
		return phenotype_method_ls
	
	@jsonify
	def getPhenotypeMethodLsJson(self):
		"""
		2009-4-6
			return a list of phenotype ids/names for html select widget 
		"""
		result = {
				'options': [
						dict(id=value, value=id) for id, value in self.getPhenotypeMethodLs(request.params.getone('call_method_id'))
						]
				}
		#result['options'].append({'id': u'[At the end]', 'value': u''})
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return result
	
	def getPhenotypeHistogramData(self, id=None):
		"""
		2009-4-6
			return histogram data for google barchart
		"""
		phenotype_method_id = request.params.get('phenotype_method_id', id)
		dc = self.getPhenotypeValue(id=phenotype_method_id, returnJson=False)
		ecotype_id2phenotype_value = dc.get("ecotype_id2phenotype_value", {})
		phenotype_value_ls = []
		for ecotype_id, phenotype_value in ecotype_id2phenotype_value.iteritems():
			if phenotype_value != 'NA':
				phenotype_value_ls.append(phenotype_value)
		import matplotlib
		count_ls, bins, patches = matplotlib.pyplot.hist(phenotype_value_ls, 20)
		
		column_name_type_ls = [("x_axis", ("string","Phenotype Value")), \
							("frequency",("number", "Frequency"))]
		return_ls = []
		no_of_bins = len(bins)
		if no_of_bins>=2:
			bin_size = abs(bins[0]-bins[1])
		else:
			bin_size = 0
		for i in range(len(count_ls)):
			x_axis = "%.2f-%.2f"%(bins[i], bins[i]+bin_size)
			frequency = count_ls[i]
			return_ls.append(dict(x_axis=x_axis, frequency=frequency))
		description = dict(column_name_type_ls)
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		response.headers['Content-Type'] = 'application/json'
		return json_result 
	
	#@jsonify
	def getPhenotypeValue(self, id=None, returnJson=True):
		"""
		2009-4-6
			given a phenotype method id
			return
				1. ecotype id : phenotype value
				2. min & max phenotype value	for requests to make phenotype markers on the map
		"""
		phenotype_method_id = request.params.get('phenotype_method_id', id)
		
		#2nd check if model.phenoData exists
		phenoData = getattr(model, 'phenoData', None)
		if phenoData is None:
			from variation.src.OutputPhenotype import OutputPhenotype
			phenoData = OutputPhenotype.getPhenotypeData(model.db.metadata.bind, phenotype_avg_table=model.Stock_250kDB.PhenotypeAvg.table.name,\
														phenotype_method_table=model.Stock_250kDB.PhenotypeMethod.table.name)
			model.phenoData = phenoData
		
		ecotype_id2phenotype_value = {}
		pheno_col_index = phenoData.col_id2col_index.get(int(phenotype_method_id))
		min_value = None
		max_value = None
		for ecotype_id,pheno_row_index in phenoData.row_id2row_index.iteritems():
			if pheno_row_index is not None and pheno_col_index is not None:
				phenotype_value = phenoData.data_matrix[pheno_row_index][pheno_col_index]
				ecotype_id2phenotype_value[ecotype_id] = phenotype_value
				if phenotype_value != 'NA':
					if min_value == None or phenotype_value<min_value:
						min_value = phenotype_value
					if max_value == None or phenotype_value>max_value:
						max_value = phenotype_value
		dc = dict(min_value=min_value, max_value=max_value, ecotype_id2phenotype_value=ecotype_id2phenotype_value)
		if returnJson:
			response.headers['Content-Type'] = 'application/json'
			return simplejson.dumps(dc)
		else:
			return dc
	
	def getPhenotypeIcon(self):
		"""
		2009-4-6
			
		"""
		value = request.params.getone('value')
		min_value = request.params.getone('min_value')
		max_value = request.params.getone('max_value')
		
		img_type = request.params.get('img_type', img_type)
		if img_type.find('svg')!=-1:
			response.headers['Content-type'] = 'image/svg'
		else:
			response.headers['Content-type'] = 'image/png'
		CandidateVsNonRatioPlot = model.Stock_250kDB.CandidateVsNonRatioPlot
		
		get_img_data_success = 0
		if id:
			row = CandidateVsNonRatioPlot.query.filter_by(results_id=results_id).\
				filter_by(type_id=type_id).filter_by(list_type_id=list_type_id).first()
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
	
	def trend(self, id=None): 
		"""
		2009-5-3
			
		"""
		c.trendDataURL = h.url_for(controller='Phenotype', action='getTrendData', id=None)
		c.trendAnnotationDataURL = h.url_for(controller='Phenotype', action='getTrendAnnotationData', id=None)
		c.multiPhenotypeDataURL = h.url_for(controller='Phenotype', action='getMultiPhenotypeData', id=None,\
										call_method_id=int(config['app_conf']['good_call_method_id']))
		return render('/PhenotypeTrend.html')
	
	def getTrendData(self):
		"""
		2009-5-3
		"""
		phenoData = hc.getPhenotypeData()
		
		column_name_type_ls = [("label", ("string","ID Name")), ("date", ("date", "Date")), \
							("which_phenotype", ("number","Which Phenotype")), ("phenotype_value", ("number","Phenotype Value"))]
		
		no_of_non_phenotype_cols = len(column_name_type_ls)
		phenotype_id_included = []
		for phenotype_id  in self.phenotype_method_id_list:
			i = phenoData.col_id2col_index.get(phenotype_id)
			if i is not None:
				phenotype_id_included.append(phenotype_id)
		
		description = dict(column_name_type_ls)
		return_ls = []
		for i in range(len(phenotype_id_included)):
			phenotype_id = phenotype_id_included[i]
			date=datetime.date(i+1000,1,1)
			for j in range(len(phenoData.row_id_ls)):
				row_label = phenoData.row_label_ls[j]
				row_id = phenoData.row_id_ls[j]
				label = '%s ID:%s'%(row_label, row_id)
				pheno_row_index = j
				pheno_col_index = phenoData.col_id2col_index.get(phenotype_id)
				if pheno_row_index is not None and pheno_col_index is not None:
					phenotype_value = phenoData.data_matrix[pheno_row_index][pheno_col_index]
					if phenotype_value == 'NA' or numpy.isnan(phenotype_value):
						phenotype_value = None
				else:
					phenotype_value = None	#numpy.nan can't be recognized by ToJSon()
				data = dict(date=date, label=label, which_phenotype=i, phenotype_value=phenotype_value)
				return_ls.append(data)
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		return json_result
	
	def getTrendAnnotationData(self):
		"""
		2009-5-3
		"""
	
	phenotype_method_id_list = range(1,8)+range(39,49) +[57,58,59, 80,81,82] + range(187,191)
	#phenotype_method_id_list = range(1,8)+ [57,58,59, 80,81,82] + range(187,191)
	phenotype_method_id_set = set(phenotype_method_id_list)
	def getMultiPhenotypeData(self):
		"""
		2009-5-3
		"""
		call_method_id = request.params.get('call_method_id', None)
		snpData = hc.getSNPDataGivenCallMethodID(call_method_id)
		phenoData = hc.getPhenotypeDataInSNPDataOrder(snpData)
		
		
		rows = model.Stock_250kDB.CallInfo.query.filter_by(method_id=call_method_id).all()
		call_info_id2pcs = {}
		for row in rows:
			call_info_id2pcs[row.id] = row.pc_value_ls[:10]
		
		column_name_type_ls = [("label", ("string","ID Name Phenotype")), ("date", ("date", "Date")), \
							("lon",("number", "Longitude")), ("lat",("number", "Latitude")), \
							("ecotypeid", ("number", "Ecotype ID")), ("name", ("string", "Native Name")), \
							("country", ("string", "Country")),\
							("pc1", ("number","PC1")), ("pc2", ("number", "PC2")), ("pc3", ("number","PC3")), \
							("pc4", ("number", "PC4")), ("pc5", ("number","PC5")), ("pc6", ("number", "PC6"))]
		
		no_of_non_phenotype_cols = len(column_name_type_ls)
		phenotype_id_included = []
		for phenotype_id  in self.phenotype_method_id_list:
			i = phenoData.col_id2col_index.get(phenotype_id)
			if i is not None:
				phenotype_name = phenoData.col_label_ls[i]
				column_name_type_ls.append(("phenotype %s"%phenotype_id, ("number", phenotype_name)))
				phenotype_id_included.append(phenotype_id)
		
		description = dict(column_name_type_ls)
		rows = model.db.metadata.bind.execute("select * from view_call where call_method_id=%s"%(call_method_id))
		return_ls = []
		phenoData_row_id_type = type(phenoData.row_id2row_index.keys()[0])
		for row in rows:
			pcs = call_info_id2pcs.get(row.call_info_id)
			if pcs:
				pc_value_ls = [getattr(call_info_pc, 'pc_value') for call_info_pc in pcs]
			else:
				pc_value_ls = [0]*10
			label = '%s ID:%s'%(row.nativename, row.ecotype_id)
			
			data = dict(date=datetime.date(2009,2,3), ecotypeid=row.ecotype_id, label=label, name=row.nativename, \
								lat=row.latitude, lon=row.longitude,\
								pc1=pc_value_ls[0], pc2=pc_value_ls[1], pc3=pc_value_ls[2], pc4=pc_value_ls[3], \
								pc5=pc_value_ls[4], pc6=pc_value_ls[5], \
								country=row.country)
			# construct a tuple with (phenotype_name, phenotype_value) as entry, incorporated into data in the end
			phenotype_name_value_tup_ls = []	#
			for i in range(no_of_non_phenotype_cols, len(column_name_type_ls)):
				phenotype_name = column_name_type_ls[i][0]
				phenotype_id = phenotype_id_included[i-no_of_non_phenotype_cols]
				pheno_row_index = phenoData.row_id2row_index.get(phenoData_row_id_type(row.ecotype_id))
				pheno_col_index = phenoData.col_id2col_index.get(phenotype_id)
				if pheno_row_index is not None and pheno_col_index is not None:
					phenotype_value = phenoData.data_matrix[pheno_row_index][pheno_col_index]
					if phenotype_value == 'NA' or numpy.isnan(phenotype_value):
						phenotype_value = None
				else:
					phenotype_value = None	#numpy.nan can't be recognized by ToJSon()
				phenotype_name_value_tup_ls.append((phenotype_name, phenotype_value))
			data.update(dict(phenotype_name_value_tup_ls))
			return_ls.append(data)
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		return json_result