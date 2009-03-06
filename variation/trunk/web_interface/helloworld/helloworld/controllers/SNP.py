import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render
from helloworld import model
import helloworld.lib.helpers as h

log = logging.getLogger(__name__)
import gviz_api	#2009-2-1 google visualization python api to export data into json-related formats
import numpy, datetime
SnpsContext = model.Stock_250kDB.SnpsContext
SNPAnnotation = model.Stock_250kDB.SNPAnnotation
Snps = model.Stock_250kDB.Snps
Gene = model.GenomeDB.Gene
Results = model.Stock_250kDB.Results
AnalysisMethod = model.Stock_250kDB.AnalysisMethod
PhenotypeMethod = model.Stock_250kDB.PhenotypeMethod
ResultsMethod = model.Stock_250kDB.ResultsMethod

class SnpController(BaseController):

	def index(self, chromosome=None, position=None, call_method_id=None, phenotype_method_id=None, analysis_method_id=None, score=-1):
		"""
		2009-2-19
			snp context
			snp annotation
			snp frequency
			geographic distribution colored by phenotype
			which accession has which allele
			phenotype distribution stratified according to allele
			significant associations (top 1k) in other phenotypes
			
		"""
		# Return a rendered template
		#   return render('/template.mako')
		# or, Return a response
		c.chromosome = request.params.get('chromosome', chromosome)
		c.position = request.params.get('position', position)
		c.call_method_id = request.params.get('call_method_id', call_method_id)
		c.phenotype_method_id = request.params.get('phenotype_method_id', phenotype_method_id)
		c.analysis_method_id = request.params.get('analysis_method_id', analysis_method_id)
		if not c.call_method_id or not c.phenotype_method_id or not c.analysis_method_id:
			c.results_id = request.params.get('results_id', None)
			if c.results_id:
				rm = ResultsMethod.get(c.results_id)
				c.call_method_id = rm.call_method_id
				c.phenotype_method_id = rm.phenotype_method_id
				c.analysis_method_id = rm.analysis_method_id
		c.score = request.params.get('score', score)
		return render('/snp.html')
	
	def getSNPSummaryInfo(self, id=None):
		"""
		2009-3-2
			information for the SNP summary table
		"""
		chromosome = request.params.get('chromosome', None)
		position = request.params.get('position', None)
		call_method_id = request.params.get('call_method_id', None)
		phenotype_method_id = request.params.get('phenotype_method_id', None)
		analysis_method_id = request.params.get('analysis_method_id', None)
		score = request.params.get('score', -1)
		
		#3rd finally construct the full data and turn it into json
		column_name_type_ls = [("chromosome", ("string", "Chromosome")), ("position", ("string", "Position")), \
							("call_method_id", ("string", "Call Method ID")), ("phenotype", ("string","Phenotype")), \
							("analysis",("string", "Association Method")), ("score", ("number","Score")), ("type_of_snp",("string", "Type of SNP")), \
							("gene_on_the_left", ("number","Gene on the Left")), ("gene_touched", ("string", "Gene touched")), ("gene_on_the_right", ("number", "Gene on the Right"))]
		
		description = dict(column_name_type_ls)
		
		c.gene_desc_names = ['gene_symbol', 'description', 'type_of_gene', 'dbxrefs']
		
		snp = Snps.query.filter_by(chromosome=chromosome).filter_by(position=position).filter(Snps.end_position==None).first()
		#raise
		#snps_context = SnpsContext.query.filter_by(snps_id=snp.snps_id).filter_by(gene_id=snp.gene_id).first()
		#row.left_or_right = getattr(snps_context, 'left_or_right', '')
		#row.disp_pos_comment = getattr(snps_context, 'disp_pos_comment', '')
		
		snp_annotation_text_ls = []
		snp_annotation_ls = SNPAnnotation.query.filter_by(snps_id=snp.id).all()
		gene_touched = None
		for snp_annotation in snp_annotation_ls:
			snp_annotation_text = snp_annotation.snp_annotation_type.short_name
			if snp_annotation.gene_id:
				snp_annotation_text += ':gene %s'%snp_annotation.gene_id
				gene = Gene.get(snp_annotation.gene_id)
				#gene_touched = '%s %s'%(gene.gene_id, gene.gene_symbol)
				gene_touched = '<a href='+h.NCBIGeneDBURL%gene.gene_id+'>%s %s</a>'%\
					(gene.gene_id, gene.gene_symbol)
			if snp_annotation.comment:
				snp_annotation_text += ':%s'%snp_annotation.comment
			if len(snp_annotation_text_ls)>0:
				if snp_annotation_text!=snp_annotation_text_ls[-1]:
					snp_annotation_text_ls.append(snp_annotation_text)
			else:
				snp_annotation_text_ls.append(snp_annotation_text)
		type_of_snp = '; '.join(snp_annotation_text_ls)
		
		phenotype_method = PhenotypeMethod.get(phenotype_method_id)
		analysis_method = AnalysisMethod.get(analysis_method_id)
		phenotype = '%s %s'%(phenotype_method.id, phenotype_method.short_name)
		
		analysis = '%s %s'%(analysis_method.id, analysis_method.short_name)
		
		return_ls = [dict(chromosome=chromosome, position=position, call_method_id=call_method_id, \
								phenotype=phenotype,\
								analysis=analysis, score=float(score), type_of_snp=type_of_snp, \
								gene_on_the_left=None, gene_touched=gene_touched, gene_on_the_right=None)]
		
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		return json_result
	
	def getSignificantHits(self, id=None):
		"""
		2009-3-2
			significant hits (stored in Stock_250kDB.Results, assume top 1000) in all phenotypes for the given SNP
			and feed back to the javascript request
		"""
		chromosome = request.params.get('chromosome', None)
		position = int(request.params.get('position', 0))
		call_method_id = request.params.get('call_method_id', None)
		phenotype_method_id = request.params.get('phenotype_method_id', None)
		analysis_method_id = request.params.get('analysis_method_id', None)
		score = request.params.get('score', None)
		
		column_name_type_ls = [("results_id", ("number", "Result Id")), ("phenotype_method_id", ("number", "Phenotype ID")), ("phenotype_short_name", ("string", "Phenotype Name")), \
							("analysis", ("string", "Association Method")), \
							("gbrowseLink", ("string","GBrowse")), ("score", ("number","Score")), \
							("rank", ("number","Rank")), ('beta', ('number', 'Beta')),\
							("maf", ("number","MAF")), ("mac", ("number","MAC")),\
							('genotype_var_perc', ('number', 'variance explained'))]
		
		description = dict(column_name_type_ls)
				
		snp = Snps.query.filter_by(chromosome=chromosome).filter_by(position=position).filter(Snps.end_position==None).first()
		
		rows = Results.query.filter_by(snps_id=snp.id).all()
		return_ls =[]
		start_pos = position-10000
		stop_pos = position+10000
		for row in rows:
			phenotype_method = row.result.phenotype_method
			analysis_method = row.result.analysis_method
			
			analysis = '%s %s'%(analysis_method.id, analysis_method.short_name)
			track_id = '%s_%s_%s'%(row.result.call_method.id, phenotype_method.id, analysis_method.id)
			gbrowseLink = "<a href=%s>%s</a>"%(h.GBrowseURL%(start_pos, stop_pos, chromosome)+track_id+"-"+track_id+"_SNP", row.result.id)
			
			return_ls.append(dict(results_id=row.results_id, phenotype_method_id=phenotype_method.id, \
								phenotype_short_name=phenotype_method.short_name,\
								analysis=analysis, gbrowseLink=gbrowseLink, \
								score=row.score, rank=row.rank, \
								beta=row.beta, maf=row.maf,\
								mac=row.mac, genotype_var_perc=row.genotype_var_perc))
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		#column_ls.sort()	#['date', 'ecotypeid', 'label', 'lat', 'lon', 'name', 'pc1', 'pc2', 'phenotype']
		json_result = data_table.ToJSon(columns_order=column_ls)
		return json_result