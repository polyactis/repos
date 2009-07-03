import logging

from helloworld.lib.base import *

log = logging.getLogger(__name__)

from formencode import htmlfill
import simplejson, re
from DisplayTopSNPTestRM import DisplaytopsnptestrmController
EntrezgeneMapping = model.GenomeDB.EntrezgeneMapping
Gene = model.GenomeDB.Gene
GeneList = model.Stock_250kDB.GeneList
GeneListType = model.Stock_250kDB.GeneListType
Results = model.Stock_250kDB.Results
ResultsMethod = model.Stock_250kDB.ResultsMethod
ResultsGene = model.Stock_250kDB.ResultsGene
ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
SnpsContext = model.Stock_250kDB.SnpsContext
SNPAnnotation = model.Stock_250kDB.SNPAnnotation
Snps = model.Stock_250kDB.Snps
from HelpOtherControllers import HelpothercontrollersController as hc
import gviz_api	#2009-2-1 google visualization python api to export data into json-related formats
import numpy, datetime, urllib
from sqlalchemy import asc

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
		2009-5-26
			handle the situation when query to table ResultsMethod returns None based on call method, phenotype, analysis.
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
		type_id = request.params.get('type_id', type_id)
		list_type_id = int(request.params.get('list_type_id', list_type_id))
		max_rank = request.params.get('max_rank', max_rank)
		if max_rank is not None:
			max_rank = int(max_rank)
		c.max_rank = max_rank
		c.type = ScoreRankHistogramType.get(type_id)
		c.list_type = model.Stock_250kDB.GeneListType.get(list_type_id)
		
		c.snp_gene_association_id2desc = self.snp_gene_association_id2desc
		c.row_ls = []
		c.counter = 0
		c.gene_desc_names = ['gene_symbol', 'description', 'type_of_gene', 'dbxrefs']
				
		if not results_id:
			call_method_id = request.params.get('call_method_id', None)
			phenotype_method_id = request.params.get('phenotype_method_id', None)
			analysis_method_id = request.params.get('analysis_method_id', None)
			rm = ResultsMethod.query.filter_by(call_method_id=call_method_id).filter_by(phenotype_method_id=phenotype_method_id).\
					filter_by(analysis_method_id=analysis_method_id).first()
			if rm is None:
				c.result = None
				return render('/display_results_gene_one.html')
			results_id = rm.id

		if results_id is None:
			results_id = request.params.get('results_id', None)
		if results_id is None:
			return 'Nothing'
		
		c.result = ResultsMethod.get(results_id)

		
		c.fetchResultsGeneURL = h.url_for(controller="DisplayResultsGene", action='fetchTopCandidateGenesFromOneResultOneGeneList', \
										results_id=results_id, type_id=type_id, list_type_id=list_type_id, max_rank=max_rank)
		
		"""
		from variation.src.GeneListRankTest import GeneListRankTest
		if list_type_id>0:	#2009-2-22
			candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
		else:
			candidate_gene_set = set()
		rows = ResultsGene.query.filter_by(results_id=results_id).\
			filter(ResultsGene.types.any(id=type_id)).\
			filter(ResultsGene.rank<=max_rank).order_by(ResultsGene.rank).order_by(ResultsGene.snps_id)
		
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
		"""
		return render('/display_results_gene_one.html')
	
	def fetchTopCandidateGenesFromOneResultOneGeneList(self, id=None, type_id=None, list_type_id=0, \
			max_rank=200):
		"""
		2009-7-2 split out of showTopCandidateGenesFromOneResultOneGeneList()
			a server API which returns a google visualization table in json
			
			used by templates underlying both showTopCandidateGenesFromOneResultOneGeneList and showResultsGeneForOnePhenotype
		"""
		results_id = request.params.get('id', id)
		type_id = request.params.get('type_id', type_id)
		list_type_id = int(request.params.get('list_type_id', list_type_id))
		max_rank = request.params.get('max_rank', max_rank)
		if max_rank is not None:
			max_rank = int(max_rank)
		c.max_rank = max_rank
		c.type = ScoreRankHistogramType.get(type_id)
		c.list_type = model.Stock_250kDB.GeneListType.get(list_type_id)
		
		c.snp_gene_association_id2desc = self.snp_gene_association_id2desc
		c.row_ls = []
		c.counter = 0
		c.gene_desc_names = ['gene_symbol', 'description', 'type_of_gene', 'dbxrefs']
		
		if not results_id:
			call_method_id = request.params.get('call_method_id', None)
			if call_method_id is None:
				call_method_id = c.type.call_method_id
			phenotype_method_id = request.params.get('phenotype_method_id', None)
			analysis_method_id = request.params.get('analysis_method_id', None)
			rm = ResultsMethod.query.filter_by(call_method_id=call_method_id).filter_by(phenotype_method_id=phenotype_method_id).\
					filter_by(analysis_method_id=analysis_method_id).first()
			if rm is None:
				c.result = None
				return "No association result found"
			results_id = rm.id

		if results_id is None:
			results_id = request.params.get('results_id', None)
		if results_id is None:
			return 'Nothing'
		
		c.result = ResultsMethod.get(results_id)
		
		from variation.src.GeneListRankTest import GeneListRankTest
		if list_type_id>0:	#2009-2-22
			candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
		else:
			candidate_gene_set = None
		return self.getAssociationsGivenGene(type_id, maxRank=c.max_rank, results_id=results_id, candidate_gene_set=candidate_gene_set)
	
	@staticmethod
	def getPhenotypeMethodLsGivenType(type_id, extra_table_name=None):
		"""
		2009-4-13
			added results_gene2type into the condition to restrict phenotype_methods
			type_id is actually used now.
		2009-3-4
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name
		type = ScoreRankHistogramType.get(type_id)
		call_method_id = type.call_method_id
		if not extra_table_name:
			extra_table_name = model.Stock_250kDB.ResultsGene.table.name
		extra_tables = ' %s c, %s y'%(extra_table_name, "results_gene2type")
		extra_condition = 'c.results_id=s.id and s.call_method_id=%s and y.score_rank_histogram_type_id=%s and y.results_gene_id=c.id'%\
				(call_method_id, type_id)
		phenotype_info  = hc.getPhenotypeInfo(affiliated_table_name=affiliated_table_name, extra_condition=extra_condition,
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
		2009-4-13
			added results_gene2type into the condition to restrict analysis_methods
			type_id is actually used now.
		2009-3-4
		"""
		affiliated_table_name = model.Stock_250kDB.ResultsMethod.table.name	#alias is 's'
		type = ScoreRankHistogramType.get(type_id)
		call_method_id = type.call_method_id
		if not extra_table_name:
			extra_table_name = model.Stock_250kDB.ResultsGene.table.name
		extra_tables = ' %s c, %s y '%(extra_table_name, "results_gene2type")
		extra_condition = 'c.results_id=s.id and s.call_method_id=%s and s.phenotype_method_id=%s and y.score_rank_histogram_type_id=%s and y.results_gene_id=c.id'%\
			(call_method_id, phenotype_method_id, type_id)
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
		list_info = hc.getListTypeInfo(affiliated_table_name=None, extra_condition=None, extra_tables=None)
		
		ls = [[0, u'No candidate gene list. (All SNPs)']]
		for i in range(len(list_info.list_type_id_ls)):
			id = list_info.list_type_id_ls[i]
			label = list_info.list_type_label_ls[i]
			ls.append([id, label])
		return ls
	
	#@jsonify
	@classmethod
	def getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethodJson(cls):
		"""
		2009-3-4
			called in the javascript-capable form
		"""
		type_id = request.params.get('type_id')
		phenotype_method_id = request.params.get('phenotype_method_id')
		analysis_method_id = request.params.get('analysis_method_id')
		if type_id and phenotype_method_id and analysis_method_id:
			type = ScoreRankHistogramType.get(type_id)
			call_method_id = type.call_method_id
			rm = model.Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).\
				filter_by(phenotype_method_id=phenotype_method_id).filter_by(analysis_method_id=analysis_method_id).first()
			results_id = rm.id
		else:
			results_id = -1
		result = {
				'results_id': results_id,
				'options': [
						dict(id=value, value=id) for id, value in cls.getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethod()
						]
				}
		#result['options'].insert(0, {'id': u'No candidate gene list. (All SNPs)', 'value': 0})
		#result['options'].insert(0, {'id': u'Please Choose ...', 'value': -1})
		return simplejson.dumps(result)
	
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
	
	def snp_gene_association_Json(self):
		"""
		2009-6-2
			a json wrap up around snp_gene_association_id_desc_ls()
		"""
		result = {
				'options': [
						dict(id=value, value=id) for id, value in self.snp_gene_association_id_desc_ls
						]
				}
		#result['options'].append({'id': u'[At the end]', 'value': u''})
		result['options'].insert(0, {'id': u'Please Choose ...', 'value': 0})
		response.headers['Content-Type'] = 'application/json'
		return simplejson.dumps(result)
	
	
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
		c.max_rank=50
		html = render('/display_results_gene_form.html')
		return htmlfill.render(html, defaults)
	
	def findGeneNameLike(self, namelike):
		"""
		2009-6-12
			limit the query return to <= 100 entries, ordered by gene symbol
		2009-4-3
			find all accessions with name beginned with 'namelike'
		"""
		tax_id = int(config['app_conf']['tax_id'])
		Gene_symbol2id = model.GenomeDB.Gene_symbol2id
		rows = Gene_symbol2id.query.filter_by(tax_id=tax_id).filter(Gene_symbol2id.table.c.gene_symbol.op('regexp')('^%s'%namelike)).\
			order_by(asc(Gene_symbol2id.gene_symbol)).limit(100)
		name_set = set()
		for row in rows:
			name_set.add(row.gene_symbol)
		return list(name_set)
	
	def geneNameAutoComplete(self):
		"""
		2009-6-2
			parameter namelike has to be of >1 characters long 
		2009-5-26
			autoComplete given a part of a gene name
		"""
		namelike = request.params.get('namelike')
		if len(namelike)<2:
			name_ls = []
		else:
			name_ls = self.findGeneNameLike(namelike)
			name_ls.sort()
			if len(name_ls)>100:
				name_ls = name_ls[:100]
		response.headers['Content-Type'] = 'application/json'
		return simplejson.dumps(dict(result=name_ls))
	
	def geneForm(self, id=None):
		"""
		2009-5-26
			form to allow search GWAS results by gene name
		"""
		if config['app_conf']['site_public']=='false' or config['app_conf']['site_public']==False:
			c.site_public = 0
		else:
			c.site_public = 1
		c.snp_gene_association_type_id = int(config['app_conf']['snp_gene_association_type_id']);
		c.snpGeneAssociationTypeURL = h.url_for(controller="DisplayResultsGene", action="snp_gene_association_Json");
		c.geneNameAutoCompleteURL = h.url_for(controller="DisplayResultsGene", action="geneNameAutoComplete");
		c.searchURL = h.url_for(controller="DisplayResultsGene", action="searchByGeneName");
		return render('/SearchGWAS.html')
	
	def findGeneIDByName(self, id=None, geneName=None):
		"""
		2009-5-26
			based on geneName, find the corresponding gene ID
		"""
		tax_id = int(config['app_conf']['tax_id'])
		Gene_symbol2id = model.GenomeDB.Gene_symbol2id
		rows = Gene_symbol2id.query.filter_by(tax_id=tax_id).filter_by(gene_symbol=geneName)
		gene_id_set = set()
		gene_id_symbol_set = set()
		for row in rows:
			if row.symbol_type=='symbol':
				gene_id_symbol_set.add(row.gene_id)
			gene_id_set.add(row.gene_id)
		if len(gene_id_set)>1:
			if len(gene_id_symbol_set)==1:
				return list(gene_id_symbol_set)[0]
			elif len(gene_id_symbol_set)==0:
				return list(gene_id_set)[0]
			elif len(gene_id_symbol_set)>1:
				return list(gene_id_symbol_set)[0]
		elif len(gene_id_set)==1:
			return list(gene_id_set)[0]
		elif len(gene_id_set)==0:
			return None
	
	def getGeneInformation(self, gene_id):
		"""
		2009-6-12
			add 'chromosome', "strand", "start", "stop" to the returned information
		2009-5-27
			return google data table in json structure of information of a gene
		"""
		column_name_type_ls = [("gene_symbol", ("string", "Symbol")), ("description", ("string", "Description")), \
							("type_of_gene", ("string", "Type Of Gene")), \
							("dbxrefs", ("string", "DB Cross References")), \
							("gene_id", ("string","Gene ID")), ("chromosome", ("string","chromosome")),\
							("strand", ("string","Strand")),\
							("start", ("number","Start")),("stop", ("number","Stop"))]
		description = dict(column_name_type_ls)
		Gene = model.GenomeDB.Gene
		gene = Gene.get(gene_id)
		EntrezgeneMapping = model.GenomeDB.EntrezgeneMapping
		em = EntrezgeneMapping.get(gene_id)
		entry = {}
		if gene is not None:
			for column_name_type in column_name_type_ls:
				column_name = column_name_type[0]
				column_type = column_name_type[1][0]
				
				if column_type=='string':
					default_value = ''
				elif column_type =='number':
					default_value = -1
				elif column_type=='date':
					default_value = datetime.date(2050, 1, 1)
				else:
					default_value = None
				if column_name in ['chromosome', "strand", "start", "stop"]:
					obj = em
				else:
					obj = gene
				
				if column_name=='gene_id':
					column_value = getattr(obj, column_name, default_value)
					if column_value:
						column_value = "<a href=%s target='_blank'>%s</a>"%(h.NCBIGeneDBURL%gene_id, column_value)
				else:
					column_value = getattr(obj, column_name, default_value)
				entry[column_name] = column_value
		return_ls = [entry]
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		response.headers['Content-Type'] = 'application/json'
		return json_result
	
	def getAssociationsGivenGene(self, typeID, minScore=None, maxRank=None, results_id=None, gene_id=None, candidate_gene_set=None):
		"""
		2009-7-2
			add argument results_id, candidate_gene_set
			make several into keyword-based ones
		2009-5-27
			return all associations related to a gene in table results_gene
		"""
		#construct the full data and turn it into json
		column_name_type_ls1 = [("results_id", ("number", "Result Id")), ("phenotype_method_id", ("string", "Phenotype ID")), \
							("phenotype_short_name", ("string", "Phenotype Name")), \
							("analysis", ("string", "Association Method"))] 
		column_name_type_ls2 = [("snps_id", ("string","SNP ID")), ("chromosome", ("string","Chr")), \
							("position", ("string","Position")), ("left_or_right", ("string","Left/Right")), \
							("disp_pos", ("string","Distance")), ("disp_pos_comment", ("string","Distance comment")), \
							("snp_annotation", ("string","SNP Annotation")), \
							("score", ("number","Score")), \
							("rank", ("number","Rank")), ('beta', ('number', 'Beta')),\
							("maf", ("number","MAF")), ("mac", ("number","MAC")),\
							('genotype_var_perc', ('number', 'variance explained'))]
		
		if gene_id is not None:
			column_name_type_ls = column_name_type_ls1 + column_name_type_ls2
			results_gene_related_columns = column_name_type_ls
		else:
			gene_column_name_type_ls = [("gene_symbol", ("string", "Symbol")), ("description", ("string", "Description")), \
								("type_of_gene", ("string", "Type Of Gene")), \
								("dbxrefs", ("string", "DB Cross References")), \
								("gene_id", ("string","Gene ID")), \
								("strand", ("string","Strand")),\
								("start", ("number","Start")),("stop", ("number","Stop"))]			
			column_name_type_ls = gene_column_name_type_ls + column_name_type_ls2
			results_gene_related_columns = column_name_type_ls2

		description = dict(column_name_type_ls)
		
		
		query = ResultsGene.query.filter(ResultsGene.types.any(id=typeID))
		if maxRank:
			query = query.filter(ResultsGene.rank<=maxRank)
		if minScore:
			query = query.filter(ResultsGene.score>=minScore)
		if results_id:
			query = query.filter_by(results_id=results_id)
		if gene_id:
			query = query.filter_by(gene_id=gene_id)
		query = query.order_by(ResultsGene.rank).order_by(ResultsGene.snps_id)
		
		return_ls = []
		counter = 0
		for row in query:
			entry = dict()
			snps_context = None
			snp_based_association_result = None
			if candidate_gene_set is not None and row.gene_id not in candidate_gene_set:	#empty set means no filtering
				continue
			if column_name_type_ls!=results_gene_related_columns:
				gene = Gene.get(row.gene_id)
				em = EntrezgeneMapping.get(row.gene_id)
				if gene is not None:
					for column_name_type in column_name_type_ls:
						column_name = column_name_type[0]
						column_type = column_name_type[1][0]
						
						if column_type=='string':
							default_value = ''
						elif column_type =='number':
							default_value = -1
						elif column_type=='date':
							default_value = datetime.date(2050, 1, 1)
						else:
							default_value = None
						if column_name in ['chromosome', "strand", "start", "stop"]:
							obj = em
						else:
							obj = gene
						
						if column_name=='gene_id':
							column_value = getattr(obj, column_name, default_value)
							if column_value:
								column_value = "<a href=%s target='_blank'>%s</a>"%(h.NCBIGeneDBURL%gene_id, column_value)
						else:
							column_value = getattr(obj, column_name, default_value)
						entry[column_name] = column_value
			
			for column_name_type in results_gene_related_columns:
				column_name = column_name_type[0]
				column_type = column_name_type[1][0]
				
				if column_type=='string':
					default_value = ''
				elif column_type =='number':
					default_value = -1
				elif column_type=='date':
					default_value = datetime.date(2050, 1, 1)
				else:
					default_value = None
				
				if column_name == 'phenotype_method_id':
					column_value = "<a href=%s target='_blank'>%s</a>"%\
						(h.url_for(controller="DisplayResults", action="showGWA", phenotype_method_id=row.result.phenotype_method_id, call_method_id=row.result.call_method_id),\
						row.result.phenotype_method_id)
				elif column_name == 'phenotype_short_name':
					column_value = row.result.phenotype_method.short_name
				elif column_name == 'analysis':
					column_value = row.result.analysis_method.short_name
				elif column_name == 'left_or_right':
					if snps_context is None:
						snps_context = SnpsContext.query.filter_by(snps_id=row.snps_id).filter_by(gene_id=row.gene_id).first()
					column_value = getattr(snps_context, column_name, default_value)
				elif column_name == 'disp_pos_comment':
					if snps_context is None:
						snps_context = SnpsContext.query.filter_by(snps_id=row.snps_id).filter_by(gene_id=row.gene_id).first()
					column_value = getattr(snps_context, column_name, default_value)
				elif column_name == 'chromosome':
					column_value = row.snp.chromosome
				elif column_name == 'position':
					column_value = row.snp.position
				elif column_name in ['beta', "maf", "mac", 'genotype_var_perc']:
					if snp_based_association_result is None:
						snp_based_association_result = Results.query.filter_by(snps_id=row.snps_id).filter_by(results_id=row.results_id).first() 
					column_value = getattr(snp_based_association_result, column_name, default_value)
				elif column_name == 'snp_annotation':
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
					column_value = ';'.join(snp_annotation_text_ls)
				elif column_name == 'snps_id':
					SNPURL = h.url_for(controller='SNP', action=None, phenotype_method_id=row.result.phenotype_method.id, \
									call_method_id=row.result.call_method.id, analysis_method_id=row.result.analysis_method.id,\
									chromosome=row.snp.chromosome, position=row.snp.position, score=row.score)
					column_value = "<a href=%s target='_blank'>%s</a>"%(SNPURL, row.snps_id)
				else:
					column_value = getattr(row, column_name, default_value)
				
				entry[column_name] = column_value
				counter += 1
			return_ls.append(entry)
		
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		response.headers['Content-Type'] = 'application/json'
		return json_result
	
	all_number_p = re.compile(r'^\d+$')
	def searchByGeneName(self, id=None):
		"""
		2009-5-26
			get association results based on gene name
		"""
		typeID = request.params.get('typeID')
		if typeID:
			typeID = int(typeID)
		minScore = request.params.get("minScore")
		if minScore:
			minScore = float(minScore)
		maxRank = request.params.get("maxRank")
		if maxRank:
			maxRank = float(maxRank)
		geneName = request.params.get('geneName')
		geneName = urllib.unquote(geneName)	#Replace %xx escapes by their single-character equivalent.
		if self.all_number_p.match(geneName):
			gene_id = int(geneName)
		else:
			gene_id = self.findGeneIDByName(geneName=geneName)
		
		mode = request.params.get('mode')
		if mode=="1":
			return self.getGeneInformation(gene_id)
		else:
			return self.getAssociationsGivenGene(typeID, minScore, maxRank, gene_id=gene_id)
	
	def showResultsGeneForOnePhenotype(self, id=None):
		"""
		2009-7-2
			similar to DisplayResults.showGWA() but with tables of genes & snps instead
		"""
		c.phenotype_method_id = int(request.params.get('phenotype_method_id', 1))
		c.call_method_id = int(request.params.get('call_method_id', config['app_conf']['published_call_method_id']))
		c.phenotypeSummaryURL = h.url_for(controller="Phenotype", action=None, phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		c.GWABaseURL = h.url_for(controller='DisplayResults', action='fetchOne', phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		c.SNPBaseURL = h.url_for(controller='SNP', action=None, phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		c.getAnalysisMethodLsURL = h.url_for(controller='DisplayResults', action='getAnalysisMethodLsJson', phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		pm = model.Stock_250kDB.PhenotypeMethod.get(c.phenotype_method_id)
		c.phenotype_method_short_name = pm.short_name
		c.phenotype_method_description = pm.method_description
		
		c.callInfoURL = h.url_for(controller='DisplayResults', action='fetchCallInfoData', id=None,\
							phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		c.phenotypeHistImageURL = h.url_for(controller='DisplayResults', action='getPhenotypeHistImage', id=None, \
										phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		c.callPhenotypeQQImageURL = h.url_for(controller='DisplayResults', action='getCallPhenotypeQQImage', id=None,\
											phenotype_method_id=c.phenotype_method_id, call_method_id=c.call_method_id)
		c.phenotypeHistogramDataURL = h.url_for(controller='Phenotype', action='getPhenotypeHistogramData', id=c.phenotype_method_id)
		
		if config['app_conf']['site_public']=='false' or config['app_conf']['site_public']==False:
			c.site_public = 0
		else:
			c.site_public = 1
		c.snpGeneAssociationTypeID = config['app_conf']['snp_gene_association_type_id']
		
		c.snpGeneAssociationTypeLsURL = h.url_for(controller="DisplayResultsGene", action="snp_gene_association_Json")
		c.snpGeneAssociationOnChangeURL = h.url_for(controller="DisplayResultsGene", action="getPhenotypeMethodLsGivenTypeJson", id=None)
		c.phenotypeMethodOnChangeURL = h.url_for(controller="DisplayResultsGene", action="getAnalysisMethodLsGivenTypeAndPhenotypeMethodJson", id=None)
		c.geneListLsURL = h.url_for(controller="DisplayResultsGene", action="getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethodJson", id=None)
		c.fetchResultsGeneURL = h.url_for(controller="DisplayResultsGene", action='fetchTopCandidateGenesFromOneResultOneGeneList', id=None)
		c.candidateGeneListURL = h.url_for(controller="DisplayResultsGene", action='showGeneList', id=None)
		
		return render('/OnePhenotypeGWASGene.html')
	
	def showGeneList(self, id=None):
		"""
		2009-7-2
			return a google visualization data table of genes based on gene list
		"""
		c.list_type_id = int(request.params.get('list_type_id', id))
		geneListType = GeneListType.get(c.list_type_id)
		c.list_short_name = geneListType.short_name
		c.list_description = geneListType.description
		c.fetchGeneListURL = h.url_for(controller="DisplayResultsGene", action='fetchGeneList', id=c.list_type_id)
		return render('/GeneListView.html')
	
	def fetchGeneList(self, id=None):
		"""
		2009-7-2
			return a google visualization data table of genes based on gene list
		"""
		list_type_id = int(request.params.get('list_type_id', id))
		column_name_type_ls = [("results_id", ("number", "Result Id")), ("phenotype_method_id", ("string", "Phenotype ID")), \
							("phenotype_short_name", ("string", "Phenotype Name")), \
							("analysis", ("string", "Association Method")), \
							("snps_id", ("string","SNP ID")), ("chromosome", ("string","Chr")), \
							("position", ("string","Position")), ("left_or_right", ("string","Left/Right")), \
							("disp_pos", ("string","Distance")), ("disp_pos_comment", ("string","Distance comment")), \
							("snp_annotation", ("string","SNP Annotation")), \
							("score", ("number","Score")), \
							("rank", ("number","Rank")), ('beta', ('number', 'Beta')),\
							("maf", ("number","MAF")), ("mac", ("number","MAC")),\
							('genotype_var_perc', ('number', 'variance explained'))]
		description = dict(column_name_type_ls)
		
		from variation.src.GeneListRankTest import GeneListRankTest
		if list_type_id>0:	#2009-2-22
			candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
		else:
			candidate_gene_set = set()
		
		query = GeneList.query.filter_by(list_type_id=list_type_id)
		if maxRank:
			query = query.filter(ResultsGene.rank<=maxRank)
		if minScore:
			query = query.filter(ResultsGene.score>=minScore)
		query = query.order_by(ResultsGene.rank).order_by(ResultsGene.snps_id)
		
		return_ls = []
		counter = 0
		for row in query:
			entry = dict()
			snps_context = None
			snp_based_association_result = None
			for column_name_type in column_name_type_ls:
				column_name = column_name_type[0]
				column_type = column_name_type[1][0]
				
				if column_type=='string':
					default_value = ''
				elif column_type =='number':
					default_value = -1
				elif column_type=='date':
					default_value = datetime.date(2050, 1, 1)
				else:
					default_value = None
				
				if column_name == 'phenotype_method_id':
					column_value = "<a href=%s target='_blank'>%s</a>"%\
						(h.url_for(controller="DisplayResults", action="showGWA", phenotype_method_id=row.result.phenotype_method_id, call_method_id=row.result.call_method_id),\
						row.result.phenotype_method_id)
				elif column_name == 'phenotype_short_name':
					column_value = row.result.phenotype_method.short_name
				elif column_name == 'analysis':
					column_value = row.result.analysis_method.short_name
				elif column_name == 'left_or_right':
					if snps_context is None:
						snps_context = SnpsContext.query.filter_by(snps_id=row.snps_id).filter_by(gene_id=row.gene_id).first()
					column_value = getattr(snps_context, column_name, default_value)
				elif column_name == 'disp_pos_comment':
					if snps_context is None:
						snps_context = SnpsContext.query.filter_by(snps_id=row.snps_id).filter_by(gene_id=row.gene_id).first()
					column_value = getattr(snps_context, column_name, default_value)
				elif column_name == 'chromosome':
					column_value = row.snp.chromosome
				elif column_name == 'position':
					column_value = row.snp.position
				elif column_name in ['beta', "maf", "mac", 'genotype_var_perc']:
					if snp_based_association_result is None:
						snp_based_association_result = Results.query.filter_by(snps_id=row.snps_id).filter_by(results_id=row.results_id).first() 
					column_value = getattr(snp_based_association_result, column_name, default_value)
				elif column_name == 'snp_annotation':
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
					column_value = ';'.join(snp_annotation_text_ls)
				elif column_name == 'snps_id':
					SNPURL = h.url_for(controller='SNP', action=None, phenotype_method_id=row.result.phenotype_method.id, \
									call_method_id=row.result.call_method.id, analysis_method_id=row.result.analysis_method.id,\
									chromosome=row.snp.chromosome, position=row.snp.position, score=row.score)
					column_value = "<a href=%s target='_blank'>%s</a>"%(SNPURL, row.snps_id)
				else:
					column_value = getattr(row, column_name, default_value)
				
				entry[column_name] = column_value
				counter += 1
			return_ls.append(entry)
		
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		response.headers['Content-Type'] = 'application/json'
		return json_result