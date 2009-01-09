import logging

from helloworld.lib.base import *

log = logging.getLogger(__name__)

from DisplayTopSNPTestRM import DisplaytopsnptestrmController

class DisplayresultsgeneController(BaseController):
	"""
	2008-10-30
		controller in charge of displaying ResultsGene
	"""
	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
		rows = ScoreRankHistogramType.query.\
			filter_by(null_distribution_type_id=1).\
			filter_by(results_type=1).\
			order_by(ScoreRankHistogramType.results_type).\
			order_by(ScoreRankHistogramType.null_distribution_type_id).\
			order_by(ScoreRankHistogramType.min_distance).\
			order_by(ScoreRankHistogramType.get_closest).\
			order_by(ScoreRankHistogramType.min_MAF).\
			order_by(ScoreRankHistogramType.allow_two_sample_overlapping).\
			all()
		c.score_rank_hist_type_ls = rows
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
												phenotype_id_ls_str = '1-7,39-59,80-82,9-13,32-38,65-74,8,60-64,75-79,158-182',\
												list_type_id_ls_str='3,6,8,28,51,64,65,68,71,76,129,24,29,30,130-139')
		
	
	def showTopCandidateGenesFromOneResultOneGeneList(self, id=None, type_id=None, list_type_id=None, max_rank=None):
		"""
		2008-11-05
			add option max_rank, which would make this function only display genes whose rank below that
		2008-10-29
		"""
		results_id = request.params.get('id', id)
		type_id = request.params.get('type_id', type_id)
		list_type_id = request.params.get('list_type_id', list_type_id)
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
		
		from variation.src.GeneListRankTest import GeneListRankTest
		#candidate_gene_list = GeneListRankTest.getGeneList(list_type_id)
		#candidate_gene_set = Set(candidate_gene_list)
		candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
		rows = ResultsGene.query.filter_by(results_id=results_id).\
			filter(ResultsGene.types.any(id=type_id)).\
			order_by(ResultsGene.rank)
		c.row_ls = []
		SnpsContext = model.Stock_250kDB.SnpsContext
		c.counter = 0
		c.gene_desc_names = ['gene_symbol', 'description', 'type_of_gene', 'dbxrefs']
		Gene = model.GenomeDB.Gene
		for row in rows:
			if row.gene_id in candidate_gene_set:
				if max_rank and row.rank>max_rank:	#2008-11-05 break given max_rank
					break
				gene = Gene.get(row.gene_id)
				for gene_desc_name in c.gene_desc_names:
					setattr(row, gene_desc_name, getattr(gene, gene_desc_name, ''))	#2008-10-29 use getattr() because gene_model could be None.
				
				snps_context = SnpsContext.query.filter_by(snps_id=row.snps_id).filter_by(gene_id=row.gene_id).first()
				row.left_or_right = getattr(snps_context, 'left_or_right', '')
				row.disp_pos_comment = getattr(snps_context, 'disp_pos_comment', '')
				
				snp_region_plot_ls = model.Stock_250kDB.SNPRegionPlot.query.\
					filter_by(center_snp_position=row.snp.position).\
					filter_by(chromosome=row.snp.chromosome).\
					filter_by(phenotype_method_id=c.result.phenotype_method_id).all()
				row.snp_region_plot_id_ls = [str(s.id) for s in snp_region_plot_ls]
				c.row_ls.append(row)
				c.counter += 1
		
		return render('/display_results_gene_one.html')