import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render, config, h, model
#from helloworld import model

from HelpOtherControllers import HelpothercontrollersController as hc
import simplejson, gviz_api, datetime

log = logging.getLogger(__name__)

class AssociationoverlapController(BaseController):

	def index(self):
		"""
		2009-11-30
		"""
		# Return a rendered template
		#   return render('/template.mako')
		# or, Return a response
		c.getCallMethodLsJsonURL = h.url_for(controller="associationoverlap", action="getCallMethodLsJson")
		c.getOverlappingDataAcrossPhenotypesURL = h.url_for(controller="associationoverlap", action="getOverlappingDataAcrossPhenotypes")
		c.getNoOfTopSNPsLsJsonURL = h.url_for(controller="associationoverlap", action="getNoOfTopSNPsLsJson")
		return render("AssociationOverlap.html")
	
	def getCallMethodLsJson(self):
		"""
		2009-11-30
			return all possible call_method_ids in table association_overlapping_stat
		"""
		affiliated_table_name = model.Stock_250kDB.AssociationOverlappingStat.table.name
		list_info  = hc.getCallMethodInfo(affiliated_table_name=affiliated_table_name)
		call_method_ls = []
		for i in range(len(list_info.id_ls)):
			id = list_info.id_ls[i]
			label = list_info.label_ls[i]
			call_method_ls.append([id, label])
			
		result = {
				'options': [
						dict(id=value, value=id) for id, value in call_method_ls
						]
				}
		#result['options'].append({'id': u'[At the end]', 'value': u''})
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return simplejson.dumps(result)
	
	def getOverlappingDataAcrossPhenotypes(self, no_of_top_snps=1000):
		"""
		2009-11-30
		"""
		call_method_id = int(request.params.get('call_method_id', config['app_conf']['published_call_method_id']))
		no_of_top_snps = int(request.params.get('no_of_top_snps', no_of_top_snps))
		phenotype_info = hc.getPhenotypeInfo(affiliated_table_name=model.Stock_250kDB.AssociationOverlappingStat.table.name, \
											extra_condition="s.call_method_id=%s"%call_method_id, extra_tables=None, \
											with_category_separator=False)
		overlapping_type_info = hc.getAssociationOverlappingTypeInfo(affiliated_table_name=model.Stock_250kDB.AssociationOverlappingStat.table.name, \
																	extra_condition="s.call_method_id=%s"%call_method_id, extra_tables=None)
		
		#construct the full data and turn it into json
		column_name_type_ls = [("label", ("string","ID Name")), ("date", ("date", "Date")), ("phenotype_method_id", ("number", "Phenotype ID")), \
							("phenotype_name", ("string", "Phenotype Name")), \
							("category", ("string","Category")), ]
		no_of_non_stat_cols = len(column_name_type_ls)
		for i in range(len(overlapping_type_info.list_type_label_ls)):
			column_name_type_ls.append(('%s'%overlapping_type_info.list_type_id_ls[i], ("number", overlapping_type_info.list_type_label_ls[i])))
		return_ls = []
		
		for i in range(len(phenotype_info.phenotype_method_id2index)):
			entry = {}
			for column_name_type in column_name_type_ls[no_of_non_stat_cols:]:	# set the default value for the overlapping stats
				column_name = column_name_type[0]
				column_type = column_name_type[1][0]
				entry[column_name] = None
			return_ls.append(entry)
		
		rows = model.Stock_250kDB.AssociationOverlappingStat.query.filter_by(call_method_id=call_method_id).filter_by(no_of_top_snps=no_of_top_snps)
		
		for row in rows:
			phenotype_index = phenotype_info.phenotype_method_id2index[row.phenotype_method_id]
			entry = return_ls[phenotype_index]
			for column_name_type in column_name_type_ls[:no_of_non_stat_cols]:	# only the info related to phenotype part
				column_name = column_name_type[0]
				column_type = column_name_type[1][0]
				
				if column_type=='string':
					default_value = ''
				elif column_type=='date':
					default_value = datetime.date(2050, 1, 1)
				else:
					default_value = None
				
				if column_name == 'phenotype_name':
					pm = model.Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
					column_value = pm.short_name
				elif column_name=='category':
					pm = model.Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
					if pm.biology_category:
						column_value = pm.biology_category.short_name
					else:
						column_value = None
				elif column_name=='label':
					pm = model.Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
					column_value = '%s %s'%(row.phenotype_method_id, pm.short_name)
				else:
					column_value = getattr(row, column_name, default_value)
					if column_value is None:
						column_value = default_value
				entry[column_name] = column_value
			# set the overlapping stat
			column_name = '%s'%row.overlapping_type_id
			entry[column_name] = row.no_of_overlapping_snps
		
		description = dict(column_name_type_ls)
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		response.headers['Content-Type'] = 'application/json'
		return json_result
	
	def getNoOfTopSNPsLsJson(self, call_method_id=32):
		"""
		2009-12-2
			return a list of no_of_top_snps given call_method_id
		"""
		call_method_id = int(request.params.get('call_method_id', call_method_id))
		
		table_name = model.Stock_250kDB.AssociationOverlappingStat.table.name
		
		rows = model.db.metadata.bind.execute("select distinct no_of_top_snps from %s"%table_name)
		
		no_of_top_snps_ls = []
		for row in rows:
			no_of_top_snps_ls.append([row.no_of_top_snps, '%s'%row.no_of_top_snps])
			
		result = {
				'options': [
						dict(id=value, value=id) for id, value in no_of_top_snps_ls
						]
				}
		#result['options'].append({'id': u'[At the end]', 'value': u''})
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return simplejson.dumps(result)