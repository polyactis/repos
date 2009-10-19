import logging

from helloworld.lib.base import *

log = logging.getLogger(__name__)
import StringIO
from pymodule import PassingData, getListOutOfStr
from formencode import htmlfill
from pylons.decorators import jsonify	#2008-12-30
from HelpOtherControllers import HelpothercontrollersController as hc


class DisplaytopsnptestrmController(BaseController):
	"""
	2008-10-26
		controller in charge of displaying results from CandidateGeneTopSNPTestRM
	"""
	def getTypes(self):
		"""
		2008-12-30
		"""
		CandidateGeneTopSNPTestRMType = model.Stock_250kDB.CandidateGeneTopSNPTestRMType

		rows = CandidateGeneTopSNPTestRMType.query.order_by(CandidateGeneTopSNPTestRMType.results_type).\
			order_by(CandidateGeneTopSNPTestRMType.test_type_id).\
			order_by(CandidateGeneTopSNPTestRMType.null_distribution_type_id).\
			order_by(CandidateGeneTopSNPTestRMType.get_closest).\
			order_by(CandidateGeneTopSNPTestRMType.min_MAF).\
			order_by(CandidateGeneTopSNPTestRMType.allow_two_sample_overlapping).all()
		return rows
	
	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		
		"""
		rows = model.db.metadata.bind.execute("select y.get_closest, y.min_MAF, y.allow_two_sample_overlapping, \
			y.results_type, y.test_type_id, y.null_distribution_type_id, group_concat(y.id order by y.id) as type_id_ls from %s y \
			group by y.get_closest, y.min_MAF, y.allow_two_sample_overlapping, \
			y.results_type, y.test_type_id, y.null_distribution_type_id"%CandidateGeneTopSNPTestRMType.table.name)
		"""
		c.rows = self.getTypes()
		"""
		for row in rows:
			row.null_distribution_type = model.Stock_250kDB.NullDistributionType.get(row.null_distribution_type_id)
			row.test_type = model.Stock_250kDB.AnalysisMethod.get(row.test_type_id)
			c.rows.append(row)
		"""
		return render('/display_top_snp_test_rm.html')
	
	@staticmethod
	def getCallMethodLsGivenType(type_id, extra_table_name=None):
		"""
		2009-2-22
			add argument extra_table_name
		2008-12-30
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name
		if not extra_table_name:
			extra_table_name = model.Stock_250kDB.CandidateGeneTopSNPTestRM.table.name
		extra_tables = ' %s c '%extra_table_name
		extra_condition = 'c.results_id=s.id and c.type_id=%s'%type_id
		list_info  = hc.getCallMethodInfo(affiliated_table_name=affiliated_table_name, extra_condition=extra_condition,
											extra_tables=extra_tables)
		call_method_ls = []
		for i in range(len(list_info.id_ls)):
			id = list_info.id_ls[i]
			label = list_info.label_ls[i]
			call_method_ls.append([id, label])
		return call_method_ls
	
	@jsonify
	def getCallMethodLsGivenTypeJson(self):
		"""
		2009-2-22
			handle one more argument: extra_table_name 
		"""
		type_id = request.params.getone('type_id')
		extra_table_name = request.params.get('extra_table_name', None)
		result = {
				'options': [
						dict(id=value, value=id) for id, value in self.getCallMethodLsGivenType(type_id, extra_table_name)
						]
				}
		#result['options'].append({'id': u'[At the end]', 'value': u''})
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return result
	
	@staticmethod
	def getPhenotypeMethodLsGivenType(type_id, call_method_id, extra_table_name=None):
		"""
		2009-2-22
			add argument extra_table_name
		2008-12-30
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name
		if not extra_table_name:
			extra_table_name = model.Stock_250kDB.CandidateGeneTopSNPTestRM.table.name
		extra_tables = ' %s c '%extra_table_name
		extra_condition = 'c.results_id=s.id and c.type_id=%s and s.call_method_id=%s'%(type_id, call_method_id)
		phenotype_info  = hc.getPhenotypeInfo(affiliated_table_name=affiliated_table_name, extra_condition=extra_condition,
											extra_tables=extra_tables)
		phenotype_method_ls = []
		for i in range(len(phenotype_info.phenotype_method_id_ls)):
			phenotype_method_id = phenotype_info.phenotype_method_id_ls[i]
			phenotype_method_label = phenotype_info.phenotype_method_label_ls[i]
			phenotype_method_ls.append([phenotype_method_id, phenotype_method_label])
		return phenotype_method_ls
	
	@jsonify
	def getPhenotypeMethodLsGivenTypeJson(self):
		"""
		2009-2-22
			
		"""
		type_id = request.params.getone('type_id')
		call_method_id = request.params.getone('call_method_id')
		extra_table_name = request.params.get('extra_table_name', None)
		result = {
				'options': [
						dict(id=value, value=id) for id, value in self.getPhenotypeMethodLsGivenType(
																									type_id, call_method_id, extra_table_name)
						]
				}
		#result['options'].append({'id': u'[At the end]', 'value': u''})
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return result
	
	@staticmethod
	def getAnalysisMethodLsGivenTypeAndPhenotypeMethod(type_id, call_method_id, phenotype_method_id, extra_table_name=None):
		"""
		2009-2-22
			add argument extra_table_name
		2008-12-30
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name	#alias is 's'
		if not extra_table_name:
			extra_table_name = model.Stock_250kDB.CandidateGeneTopSNPTestRM.table.name
		extra_tables = ' %s c '%extra_table_name
		extra_condition = 'c.results_id=s.id and c.type_id=%s and s.call_method_id=%s and s.phenotype_method_id=%s'%\
			(type_id, call_method_id, phenotype_method_id)
		list_info = hc.getAnalysisMethodInfo(affiliated_table_name, extra_condition=extra_condition, extra_tables=extra_tables)
		
		analysis_method_ls = []
		
		for i in range(len(list_info.id_ls)):
			id = list_info.id_ls[i]
			label = list_info.label_ls[i]
			analysis_method_ls.append([id, label])
		return analysis_method_ls
	
	@jsonify
	def getAnalysisMethodLsGivenTypeAndPhenotypeMethodJson(self):
		"""
		2009-2-24
			get extra_table_name from request.params
		"""
		type_id = request.params.getone('type_id')
		call_method_id = request.params.getone('call_method_id')
		phenotype_method_id = request.params.getone('phenotype_method_id')
		extra_table_name = request.params.get('extra_table_name', None)
		result = {
				'options': [
						dict(id=value, value=id) for id, value in self.getAnalysisMethodLsGivenTypeAndPhenotypeMethod(
																									type_id, call_method_id, phenotype_method_id, extra_table_name)
						]
				}
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return result
	
	@staticmethod
	def getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethod(type_id, call_method_id, phenotype_method_id, analysis_method_id, affiliated_table_name=None):
		"""
		2009-2-22
			add argument affiliated_table_name (but not used)
		2008-12-30
		"""
		if not affiliated_table_name:
			affiliated_table_name = model.Stock_250kDB.CandidateGeneTopSNPTestRM.table.name	#alias is 's'
		extra_tables = ' %s c '%model.Stock_250kDB.ResultsMethod.table.name
		extra_condition = 's.results_id=c.id and s.type_id=%s and c.call_method_id=%s and c.phenotype_method_id=%s and c.analysis_method_id=%s'%\
			(type_id, call_method_id, phenotype_method_id, analysis_method_id)
		list_info = hc.getListTypeInfo(affiliated_table_name, extra_condition=extra_condition, extra_tables=extra_tables)
		
		ls = []
		
		for i in range(len(list_info.list_type_id_ls)):
			id = list_info.list_type_id_ls[i]
			label = list_info.list_type_label_ls[i]
			ls.append([id, label])
		return ls
	
	@jsonify
	def getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethodJson(self):
		type_id = request.params.getone('type_id')
		call_method_id = request.params.getone('call_method_id')
		phenotype_method_id = request.params.getone('phenotype_method_id')
		analysis_method_id = request.params.getone('analysis_method_id')
		rm = model.Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).\
			filter_by(phenotype_method_id=phenotype_method_id).filter_by(analysis_method_id=analysis_method_id).first()
		results_id = rm.id
		result = {
				'results_id': results_id,
				'options': [
						dict(id=value, value=id) for id, value in self.getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethod(
																						type_id, call_method_id, phenotype_method_id, analysis_method_id)
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
		defaults = {'type_id': id}
		rows = self.getTypes()
		c.types = []
		type_attri_name_label_ls = [['get_closest', 'get_closest'], ['min_MAF','MAF'], ['results_type','results_type'], \
								['test_type_id', 'test type'], ['null_distribution_type_id', 'Null distribution type']]
		for row in rows:
			value_ls = []
			for attri_name, label in type_attri_name_label_ls:
				if attri_name=='test_type_id':
					attri = row.test_type.short_name
				elif attri_name=='null_distribution_type_id':
					attri = row.null_distribution_type.short_name
				else:
					attri = getattr(row, attri_name, '')
				value = '%s=%s'%(label, attri)
				value_ls.append(value)
			value = ','.join(value_ls)
			value = '%s :%s'%(row.id, value)
			c.types.append([row.id, value])
		
		c.call_method_ls = self.getCallMethodLsGivenType(id)
		c.call_method_ls.insert(0, [0, u'Please Choose ...'])
		html = render('/display_top_snp_test_rm_form.html')
		return htmlfill.render(html, defaults)
	
	@classmethod
	def type(cls, id=None, TypeClass=model.Stock_250kDB.CandidateGeneTopSNPTestRMType, \
			template='/display_top_snp_test_rm_type.html', \
			phenotype_id_ls_str = '1-7,39-61,80-82,9-13,32-38,65-74', \
			list_type_id_ls_str = '3,6,8,28,51,64,65,68,71,76,129,24,29,30,130,131',\
			analysis_method_id_ls_str = '1,5,6,7',\
			figure_type=1, call_method_id=17):
		"""
		2009-10-13
			allow user to specify phenotype_id_ls_str and list_type_id_ls_str in URL
			add figure_type argument:
				1: candidate vs non-candidate ratio (in db beforehand)
				2: 2D matrix-plot of pvalue for enrichment statistic. row is cutoff. column is distance between gene and SNP.
			add call_method_id argument
			add argument analysis_method_id_ls_str
		2008-10-30
			generalized so that DisplayResultsGene could call it
		"""
		ResultsMethod = model.Stock_250kDB.ResultsMethod
		
		if id is None:
			id = request.params.get('id', '1')
		
		phenotype_id_ls_str = request.params.get('phenotype_id_ls_str', phenotype_id_ls_str)	# 2009-10-13
		list_type_id_ls_str = request.params.get('list_type_id_ls_str', list_type_id_ls_str)	# 2009-10-13
		analysis_method_id_ls_str  = request.params.get('analysis_method_id_ls_str', analysis_method_id_ls_str)	# 2009-10-13
		c.figure_type = request.params.get('figure_type', figure_type)	# 2009-10-13
		call_method_id = request.params.get('call_method_id', call_method_id)	# 2009-10-13
		
		phenotype_id_ls = getListOutOfStr(phenotype_id_ls_str, data_type=str)
		
		list_type_id_ls = getListOutOfStr(list_type_id_ls_str, data_type=str)
		analysis_method_id_ls = getListOutOfStr(analysis_method_id_ls_str, data_type=int)
		
		#extra_tables = ' %s c '%(CandidateGeneTopSNPTestRM.table.name)
		extra_condition = 'p.id in (%s)'%(','.join(phenotype_id_ls))
		
		c.phenotype_info  = hc.getPhenotypeInfo(extra_condition=extra_condition, extra_tables=None)
		#extra_condition = 's.type_id =%s'%id
		extra_condition = 'p.id in (%s)'%(','.join(list_type_id_ls))
		
		c.list_info = hc.getListTypeInfo(extra_condition=extra_condition)
		
		data_matrix = []
		no_of_rows = len(c.list_info.list_type_id_ls)
		no_of_cols = len(c.phenotype_info.phenotype_method_id_ls)
		for i in range(no_of_rows):
			data_matrix.append([None]*no_of_cols)
		
		c.type = TypeClass.get(id)
		
		#2008-10-29 a faster way to come up the data_matrix but data is not guaranteed inside
		rows = ResultsMethod.query.filter(ResultsMethod.phenotype_method_id.in_(map(int, phenotype_id_ls))).\
			filter(ResultsMethod.analysis_method_id.in_(analysis_method_id_ls)).\
			filter_by(call_method_id=call_method_id).\
			order_by(ResultsMethod.analysis_method_id)
		counter = 0
		for row in rows:
			col_index = c.phenotype_info.phenotype_method_id2index.get(row.phenotype_method_id)
			if col_index is None:
				continue
			for list_type_id, row_index in c.list_info.list_type_id2index.iteritems():
				if list_type_id>0:
					if data_matrix[row_index][col_index] is None:
						data_matrix[row_index][col_index] = []
					data_matrix[row_index][col_index].append((row.id, row.analysis_method.short_name))
					counter +=1
		"""
		rows = CandidateGeneTopSNPTestRM.query.filter_by(type_id=id)
		
		for row in rows:
			row_index = c.list_info.list_type_id2index[row.list_type_id]
			col_index = c.phenotype_info.phenotype_method_id2index[row.result.phenotype_method_id]
			if data_matrix[row_index][col_index] is None:
				data_matrix[row_index][col_index] = Set()
			data_matrix[row_index][col_index].add(row.results_id)
			counter +=1
		"""
		c.counter = counter
		c.data_matrix = data_matrix
		return render(template)
	
	def generate2DPvalueForEnrichmentStatFigure(self, results_id=None, type_id=None, list_type_id=None, get_data_matrix=False):
		"""
		2009-10-13
			split out of showOne()
		"""
		plot_fname = 'top_snp_test_%s_%s_%s_pvalue.png'%(results_id, list_type_id, type_id)
		plot_file_path = os.path.join(config['app_conf']['plots_store'], plot_fname)
		if get_data_matrix:
			from_where_clause = "from %s t, %s y where t.type_id=y.id and t.results_id=%s and t.list_type_id=%s and y.id =%s"%\
				(model.Stock_250kDB.CandidateGeneTopSNPTestRM.table.name, model.Stock_250kDB.CandidateGeneTopSNPTestRMType.table.name,\
				results_id, list_type_id, type_id)
			from variation.src.DrawTopSNPTest2DMapForOneRM import DrawTopSNPTest2DMapForOneRM
			no_of_top_snps_info = DrawTopSNPTest2DMapForOneRM.get_no_of_top_snps_info(model.db, from_where_clause)
			min_distance_info = DrawTopSNPTest2DMapForOneRM.get_min_distance_info(model.db, from_where_clause)
			rdata = DrawTopSNPTest2DMapForOneRM.get_data_matrix(model.db, no_of_top_snps_info, min_distance_info, from_where_clause, need_other_values=True)
		else:
			no_of_top_snps_info = None
			min_distance_info = None
			rdata = None
		
		if os.path.isfile(plot_file_path):
			return_data = PassingData(no_of_top_snps_info=no_of_top_snps_info, \
									min_distance_info=min_distance_info, rdata=rdata, pvalue_png_fname=plot_fname)
			return return_data
		elif rdata is None:	# 2009-10-18 must fetch rdata in order to draw the figure.
			from_where_clause = "from %s t, %s y where t.type_id=y.id and t.results_id=%s and t.list_type_id=%s and y.id =%s"%\
				(model.Stock_250kDB.CandidateGeneTopSNPTestRM.table.name, model.Stock_250kDB.CandidateGeneTopSNPTestRMType.table.name,\
				results_id, list_type_id, type_id)
			from variation.src.DrawTopSNPTest2DMapForOneRM import DrawTopSNPTest2DMapForOneRM
			no_of_top_snps_info = DrawTopSNPTest2DMapForOneRM.get_no_of_top_snps_info(model.db, from_where_clause)
			min_distance_info = DrawTopSNPTest2DMapForOneRM.get_min_distance_info(model.db, from_where_clause)
			rdata = DrawTopSNPTest2DMapForOneRM.get_data_matrix(model.db, no_of_top_snps_info, min_distance_info, from_where_clause, need_other_values=True)

		
		from pymodule.DrawMatrix import drawMatrixLegend
		
		#draw the pvalue plot
		curve_fname = 'top_snp_test_%s_%s_%s_curve.png'%(results_id, list_type_id, type_id)
		curve_file_path = str(os.path.join(config['app_conf']['plots_store'], curve_fname))	#str() to avoid "TypeError: cannot return std::string from Unicode object"
		
		if len(rdata.data_matrix)>1 and len(rdata.data_matrix[0])>=1:
			im = drawMatrixLegend(rdata.data_matrix, left_label_ls=no_of_top_snps_info.label_ls, \
							top_label_ls=min_distance_info.label_ls, min_value=rdata.min_value, max_value=2)	# rdata.max_value
			png_data = StringIO.StringIO()
			im.save(plot_file_path, format='png')
			#DrawTopSNPTest2DMapForOneRM.plotCurve(rdata, no_of_top_snps_info, min_distance_info, curve_file_path)
			return_data = PassingData(no_of_top_snps_info=no_of_top_snps_info, \
									min_distance_info=min_distance_info, rdata=rdata, pvalue_png_fname=plot_fname)
			return return_data
		else:
			return None
	
	def getImage(self, id=None):
		"""
		2008-10-26
			get image from a file(path=id), which is generated by showOne()
		"""
		response.headers['Content-type'] = 'image/png'
		plot_file_path = os.path.join(config['app_conf']['plots_store'], id)
		img_data = None
		if os.path.isfile(plot_file_path):
			inf = open(plot_file_path, 'rb')
			img_data = inf.read()
			del inf
		if not img_data:
			response.headers['Content-type'] = 'text/plain'
			img_data = 'No Image'
		return img_data
		
	def showOne(self, id=None, type_id=None, list_type_id=None):
		"""
		2009-10-18
			some code split to generate2DPvalueForEnrichmentStatFigure(), for others to call
		2008-10-30
			add code to DrawTopSNPTest2DMapForOneRM.plotCurve()
		2008-10-24
			show all data identified by (results_id=id, type_id=type_id, list_type_id=list_type_id)
		"""
		type_id = request.params.get('type_id', type_id)
		list_type_id = request.params.get('list_type_id', list_type_id)
		results_id = request.params.get('id', id)
		if results_id is None:
			results_id = request.params.get('results_id', None)
		if results_id is None:
			return 'Nothing'
		
		c.type_id = type_id
		c.result = model.Stock_250kDB.ResultsMethod.get(results_id)
		c.list_type = model.Stock_250kDB.GeneListType.get(list_type_id)
		c.CandidateGeneTopSNPTestRMType_id_min_distance2ScoreRankHistogramType_id = model.CandidateGeneTopSNPTestRMType_id_min_distance2ScoreRankHistogramType_id
		
		return_data = self.generate2DPvalueForEnrichmentStatFigure(results_id=results_id, type_id=type_id, list_type_id=list_type_id,\
																get_data_matrix=True)
		
		if return_data:
			c.pvalue_png_fname = return_data.pvalue_png_fname
			c.no_of_top_snps_info = return_data.no_of_top_snps_info
			c.min_distance_info = return_data.min_distance_info
			c.rdata = return_data.rdata
			#c.rdata.data_matrix_non_candidate_gw_size = 250000-c.rdata.data_matrix_candidate_gw_size
		
		
		curve_fname = 'top_snp_test_%s_%s_%s_curve.png'%(results_id, list_type_id, type_id)
		curve_file_path = str(os.path.join(config['app_conf']['plots_store'], curve_fname))	#str() to avoid "TypeError: cannot return std::string from Unicode object"
		c.curve_fname = curve_fname
		
		#	self.savePlot(c.snp_region_plot, plot_file_path)
		"""
		im = drawMatrixLegend(c.rdata.data_matrix_candidate_sample_size/c.rdata.data_matrix_candidate_gw_size, \
							left_label_ls=no_of_top_snps_info.label_ls, \
							top_label_ls=min_distance_info.label_ls)
		"""
		plot_fname = 'top_snp_test_%s_%s_%s_candidate_sample.png'%(results_id, list_type_id, type_id)
		plot_file_path = os.path.join(config['app_conf']['plots_store'], plot_fname)
		#im.save(plot_file_path, format='png')
		c.candidate_sample_png_fname = plot_fname
		
		"""
		im = drawMatrixLegend(c.rdata.data_matrix_non_candidate_sample_size, left_label_ls=no_of_top_snps_info.label_ls, \
							top_label_ls=min_distance_info.label_ls)
		"""
		plot_fname = 'top_snp_test_%s_%s_%s_non_candidate_sample.png'%(results_id, list_type_id, type_id)
		plot_file_path = os.path.join(config['app_conf']['plots_store'], plot_fname)
		#im.save(plot_file_path, format='png')
		c.non_candidate_sample_png_fname = plot_fname
		
		return render('/display_top_snp_test_rm_one.html')
	
	def queryImage(self, id=None, type_id=None, list_type_id=28, img_type='png_data'):
		"""
		2009-10-13
			display plots from table CandidateVsNonRatioPlot
			img_type could be
				'png_data': png figure of candidate vs non-cand enrichment ratio
				'svg_data': svg version of above
				'pvalue2D': pvalue for each enrichment ratio
		"""
		type_id = request.params.get('type_id', type_id)
		list_type_id = request.params.get('list_type_id', list_type_id)
		results_id = request.params.get('id', id)
		
		if results_id is None:
			results_id = request.params.get('results_id', None)
		if results_id is None:
			response.headers['Content-type'] = 'text/plain'
			return 'No image'
		
		type_id = int(type_id)
		list_type_id = int(list_type_id)
		results_id = int(results_id)
		
		
		img_type = request.params.get('img_type', img_type)
		if img_type.find('svg')!=-1:
			response.headers['Content-type'] = 'image/svg'
		else:
			response.headers['Content-type'] = 'image/png'
		CandidateVsNonRatioPlot = model.Stock_250kDB.CandidateVsNonRatioPlot
		
		get_img_data_success = 0
		
		if results_id:
			if img_type=='pvalue2D':
				return_data = self.generate2DPvalueForEnrichmentStatFigure(results_id=results_id, type_id=type_id, list_type_id=list_type_id)
				if return_data:
					return self.getImage(return_data.pvalue_png_fname)
					get_img_data_success = 1
			
			else:
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
	
