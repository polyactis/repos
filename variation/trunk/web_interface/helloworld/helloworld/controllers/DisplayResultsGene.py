import logging

from helloworld.lib.base import *

log = logging.getLogger(__name__)

from formencode import htmlfill

from DisplayTopSNPTestRM import DisplaytopsnptestrmController
SnpsContext = model.Stock_250kDB.SnpsContext
SNPAnnotation = model.Stock_250kDB.SNPAnnotation
ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
ResultsMethod = model.Stock_250kDB.ResultsMethod

class DisplayresultsgeneController(BaseController):
	"""
	2008-10-30
		controller in charge of displaying ResultsGene
	"""
	def snp_gene_association_types(self):
		"""
		2009-3-4
			refactored out of index() so that form() could call it as well.
		"""
		rows = ScoreRankHistogramType.query.\
			filter_by(null_distribution_type_id=1).\
			filter_by(results_type=1).\
			order_by(ScoreRankHistogramType.call_method_id).\
			order_by(ScoreRankHistogramType.results_type).\
			order_by(ScoreRankHistogramType.null_distribution_type_id).\
			order_by(ScoreRankHistogramType.min_distance).\
			order_by(ScoreRankHistogramType.get_closest).\
			order_by(ScoreRankHistogramType.min_MAF).\
			order_by(ScoreRankHistogramType.allow_two_sample_overlapping).\
			all()
		return rows
	snp_gene_association_types = property(snp_gene_association_types)
	
	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		
		c.score_rank_hist_type_ls = self.snp_gene_association_types
		"""
		for row in rows:	#2008-10-29 only the ones with results_gene associated
			if len(row.results_gene_ls)>0:
				c.score_rank_hist_type_ls.append(row)
		"""
		return render('/display_results_gene.html')
	
	def gene_list_by_phenotype(self, id=None):
		"""
		2008-10-29
			display it as a 2-D matrix with candidate gene lists as rows and phenotypes as columns
		"""
		ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
		return DisplaytopsnptestrmController.type(id, TypeClass=ScoreRankHistogramType, \
												template='/display_results_gene_list_by_phenotype.html',\
												phenotype_id_ls_str = '1-7,39-59,80-82,9-13,14-31,32-38,65-74,8,60-64,75-79,158-186,',\
												list_type_id_ls_str='3,6,8,28,51,64,65,68,71,76,129,24,29,30,130-139,141,144,86,43')
		
	
	def showTopCandidateGenesFromOneResultOneGeneList(self, id=None, type_id=None, list_type_id=0, \
			max_rank=200):
		"""
		2009-3-4
			max_rank default is 200.
			pass snp_gene_association_id2desc to template through c
		2009-2-22
			if list_type_id (candidate_gene_list_type) is 0, no filtering, take all SNPs.
		2009-1-28
			return snp annotation to the template as well
		2008-11-05
			add option max_rank, which would make this function only display genes whose rank below that
		2008-10-29
		"""
		results_id = request.params.get('id', id)
		if not results_id:
			call_method_id = request.params.get('call_method_id', None)
			phenotype_method_id = request.params.get('phenotype_method_id', None)
			analysis_method_id = request.params.get('analysis_method_id', None)
			rm = ResultsMethod.query.filter_by(call_method_id=call_method_id).filter_by(phenotype_method_id=phenotype_method_id).\
					filter_by(analysis_method_id=analysis_method_id).first()
			results_id = rm.id
		type_id = request.params.get('type_id', type_id)
		list_type_id = int(request.params.get('list_type_id', list_type_id))
		max_rank = request.params.get('max_rank', max_rank)
		if results_id is None:
			results_id = request.params.get('results_id', None)
		if results_id is None:
			return 'Nothing'
		if max_rank is not None:
			max_rank = int(max_rank)
		
		ResultsGene = model.Stock_250kDB.ResultsGene
		ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
		c.result = model.Stock_250kDB.ResultsMethod.get(results_id)
		c.type = ScoreRankHistogramType.get(type_id)
		c.list_type = model.Stock_250kDB.GeneListType.get(list_type_id)
		
		c.snp_gene_association_id2desc = self.snp_gene_association_id2desc
		
		from variation.src.GeneListRankTest import GeneListRankTest
		if list_type_id>0:	#2009-2-22
			candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
		else:
			candidate_gene_set = set()
		rows = ResultsGene.query.filter_by(results_id=results_id).\
			filter(ResultsGene.types.any(id=type_id)).\
			filter(ResultsGene.rank<=max_rank).order_by(ResultsGene.rank).order_by(ResultsGene.snps_id)
		c.row_ls = []
		c.counter = 0
		c.gene_desc_names = ['gene_symbol', 'description', 'type_of_gene', 'dbxrefs']
		Gene = model.GenomeDB.Gene
		for row in rows:
			if row.gene_id in candidate_gene_set or list_type_id==0:	#2009-2-22
				if max_rank and row.rank>max_rank:	#2008-11-05 break given max_rank
					break
				gene = Gene.get(row.gene_id)
				for gene_desc_name in c.gene_desc_names:
					setattr(row, gene_desc_name, getattr(gene, gene_desc_name, ''))	#2008-10-29 use getattr() because gene_model could be None.
				
				snps_context = SnpsContext.query.filter_by(snps_id=row.snps_id).filter_by(gene_id=row.gene_id).first()
				row.left_or_right = getattr(snps_context, 'left_or_right', '')
				row.disp_pos_comment = getattr(snps_context, 'disp_pos_comment', '')
				
				snp_annotation_text_ls = []
				snp_annotation_ls = SNPAnnotation.query.filter_by(snps_id=row.snps_id).filter_by(gene_id=row.gene_id).all()
				for snp_annotation in snp_annotation_ls:
					snp_annotation_text = snp_annotation.snp_annotation_type.short_name
					if snp_annotation.comment:
						snp_annotation_text += ':%s'%snp_annotation.comment
					if len(snp_annotation_text_ls)>0:
						if snp_annotation_text!=snp_annotation_text_ls[-1]:
							snp_annotation_text_ls.append(snp_annotation_text)
					else:
						snp_annotation_text_ls.append(snp_annotation_text)
				row.snp_annotation = ';'.join(snp_annotation_text_ls)
				
				snp_region_plot_ls = model.Stock_250kDB.SNPRegionPlot.query.\
					filter_by(center_snp_position=row.snp.position).\
					filter_by(chromosome=row.snp.chromosome).\
					filter_by(phenotype_method_id=c.result.phenotype_method_id).all()
				row.snp_region_plot_id_ls = [str(s.id) for s in snp_region_plot_ls]
				c.row_ls.append(row)
				c.counter += 1
		c.max_rank = max_rank
		return render('/display_results_gene_one.html')
	
	@staticmethod
	def getPhenotypeMethodLsGivenType(type_id, extra_table_name=None):
		"""
		2009-3-4
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name
		type = ScoreRankHistogramType.get(type_id)
		call_method_id = type.call_method_id
		if not extra_table_name:
			extra_table_name = model.Stock_250kDB.ResultsGene.table.name
		extra_tables = ' %s c '%extra_table_name
		extra_condition = 'c.results_id=s.id and s.call_method_id=%s'%(call_method_id)
		phenotype_info  = h.getPhenotypeInfo(affiliated_table_name=affiliated_table_name, extra_condition=extra_condition,
											extra_tables=extra_tables)
		phenotype_method_ls = []
		for i in range(len(phenotype_info.phenotype_method_id_ls)):
			phenotype_method_id = phenotype_info.phenotype_method_id_ls[i]
			phenotype_method_label = phenotype_info.phenotype_method_label_ls[i]
			phenotype_method_ls.append([phenotype_method_id, phenotype_method_label])
		return phenotype_method_ls, call_method_id
	
	@jsonify
	def getPhenotypeMethodLsGivenTypeJson(self):
		"""
		2009-3-4
			called in the javascript-capable form
		"""
		type_id = request.params.get('type_id')
		extra_table_name = request.params.get('extra_table_name', None)
		phenotype_method_ls, call_method_id = self.getPhenotypeMethodLsGivenType(type_id, extra_table_name)
		result = {
				'call_method_id': call_method_id,
				'options': [
						dict(id=value, value=id) for id, value in phenotype_method_ls
						]
				}
		#result['options'].append({'id': u'[At the end]', 'value': u''})
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return result
	
	@staticmethod
	def getAnalysisMethodLsGivenTypeAndPhenotypeMethod(type_id, phenotype_method_id, extra_table_name=None):
		"""
		2009-3-4
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name	#alias is 's'
		type = ScoreRankHistogramType.get(type_id)
		call_method_id = type.call_method_id
		if not extra_table_name:
			extra_table_name = model.Stock_250kDB.ResultsGene.table.name
		extra_tables = ' %s c '%extra_table_name
		extra_condition = 'c.results_id=s.id and s.call_method_id=%s and s.phenotype_method_id=%s'%\
			(call_method_id, phenotype_method_id)
		list_info = h.getAnalysisMethodInfo(affiliated_table_name, extra_condition=extra_condition, extra_tables=extra_tables)
		
		analysis_method_ls = []
		
		for i in range(len(list_info.id_ls)):
			id = list_info.id_ls[i]
			label = list_info.label_ls[i]
			analysis_method_ls.append([id, label])
		return analysis_method_ls
	
	@jsonify
	def getAnalysisMethodLsGivenTypeAndPhenotypeMethodJson(self):
		"""
		2009-3-4
			called in the javascript-capable form
		"""
		type_id = request.params.get('type_id', 1)
		phenotype_method_id = request.params.get('phenotype_method_id')
		extra_table_name = request.params.get('extra_table_name', None)
		result = {
				'options': [
						dict(id=value, value=id) for id, value in self.getAnalysisMethodLsGivenTypeAndPhenotypeMethod(
																									type_id, phenotype_method_id, extra_table_name)
						]
				}
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		return result
	
	@staticmethod
	def getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethod():
		"""
		2009-3-4
			called by getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethodJson()
		"""
		#if not affiliated_table_name:
		#	affiliated_table_name = model.Stock_250kDB.CandidateGeneTopSNPTestRM.table.name	#alias is 's'
		#extra_tables = ' %s c '%model.Stock_250kDB.ResultsMethod.table.name
		#extra_condition = 's.results_id=c.id and s.type_id=%s and c.call_method_id=%s and c.phenotype_method_id=%s and c.analysis_method_id=%s'%\
		#	(type_id, call_method_id, phenotype_method_id, analysis_method_id)
		list_info = h.getListTypeInfo(affiliated_table_name=None, extra_condition=None, extra_tables=None)
		
		ls = [[0, u'No candidate gene list. (All SNPs)']]
		for i in range(len(list_info.list_type_id_ls)):
			id = list_info.list_type_id_ls[i]
			label = list_info.list_type_label_ls[i]
			ls.append([id, label])
		return ls
	
	@jsonify
	def getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethodJson(self):
		"""
		2009-3-4
			called in the javascript-capable form
		"""
		type_id = request.params.get('type_id')
		type = ScoreRankHistogramType.get(type_id)
		call_method_id = type.call_method_id
		phenotype_method_id = request.params.get('phenotype_method_id')
		analysis_method_id = request.params.get('analysis_method_id')
		rm = model.Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).\
			filter_by(phenotype_method_id=phenotype_method_id).filter_by(analysis_method_id=analysis_method_id).first()
		results_id = rm.id
		result = {
				'results_id': results_id,
				'options': [
						dict(id=value, value=id) for id, value in self.getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethod()
						]
				}
		#result['options'].insert(0, {'id': u'No candidate gene list. (All SNPs)', 'value': 0})
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': -1})
		return result
	
	def snp_gene_association_id_desc_ls(self):
		"""
		2009-3-4
			refactored out of form()
		"""
		snp_gene_association_id_desc_ls = []
		type_attri_name_label_ls = [['call_method_id', 'call_method_id'], ['min_distance', 'maxGeneToSNPdistance'], \
								['get_closest', 'get_closest'], ['min_MAF','MAF'], ['results_type','results_type'], \
								['allow_two_sample_overlapping', 'allow_two_sample_overlapping'], ['null_distribution_type_id', 'Null distribution type']]
		for row in self.snp_gene_association_types:
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
			value = '%s: %s'%(row.id, value)
			snp_gene_association_id_desc_ls.append([row.id, value])
		return snp_gene_association_id_desc_ls
	snp_gene_association_id_desc_ls = property(snp_gene_association_id_desc_ls)
	
	def snp_gene_association_id2desc(self):
		"""
		2009-3-4
			passed to template in showTopCandidateGenesFromOneResultOneGeneList()
		"""
		snp_gene_association_id2desc = {}
		for id, desc in self.snp_gene_association_id_desc_ls:
			snp_gene_association_id2desc[id] = desc
		return snp_gene_association_id2desc
	snp_gene_association_id2desc = property(snp_gene_association_id2desc)
	
	def form(self, id=None):
		"""
		2009-3-4
			modelled after form() in DisplayTopSNPTestRM.py
		"""
		if id is None:
			id = request.params.get('id', 1)
		defaults = {'type_id': id}
		c.types = self.snp_gene_association_id_desc_ls
		c.gene_list_ls = self.getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethod()
		c.max_rank=500
		html = render('/display_results_gene_form.html')
		return htmlfill.render(html, defaults)