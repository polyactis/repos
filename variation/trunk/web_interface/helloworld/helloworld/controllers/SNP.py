import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render
from helloworld import model
import helloworld.lib.helpers as h

log = logging.getLogger(__name__)
import gviz_api	#2009-2-1 google visualization python api to export data into json-related formats
import numpy, datetime
from variation.src.common import getEcotypeInfo
from pymodule import number2nt
from HelpOtherControllers import HelpothercontrollersController as hc


SnpsContext = model.Stock_250kDB.SnpsContext
SNPAnnotation = model.Stock_250kDB.SNPAnnotation
Snps = model.Stock_250kDB.Snps
Gene = model.GenomeDB.Gene
Results = model.Stock_250kDB.Results
AnalysisMethod = model.Stock_250kDB.AnalysisMethod
PhenotypeMethod = model.Stock_250kDB.PhenotypeMethod
ResultsMethod = model.Stock_250kDB.ResultsMethod

class SnpController(BaseController):
	def index(self, id=None):
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
		self.readInRequestParams(request, c)
		
		start_pos = c.position-40000
		stop_pos = c.position+40000
		track_id = '%s_%s_%s'%(c.call_method_id, c.phenotype_method_id, c.analysis_method_id)
		c.gbrowseLink = h.GBrowseURL%(start_pos, stop_pos, c.chromosome)+track_id+"-"+track_id+"_SNP"
		
		phenotype_method = model.Stock_250kDB.PhenotypeMethod.get(c.phenotype_method_id)
		c.pageTitle = "SNP chromosome %s position %s Phenotype %s %s"%\
			(c.chromosome, c.position, c.phenotype_method_id, phenotype_method.short_name)
		#return render('/snp.html')
		return render('/SNP.html')
	
	def readInRequestParams(self, request, param_obj):
		"""
		2009-3-5
			common function to read in parameters from http request
		"""
		param_obj.chromosome = request.params.get('chromosome', None)
		param_obj.position = request.params.get('position', None)
		param_obj.call_method_id = request.params.get('call_method_id', None)
		param_obj.phenotype_method_id = request.params.get('phenotype_method_id', None)
		param_obj.analysis_method_id = request.params.get('analysis_method_id', None)
		param_obj.score = request.params.get('score', -1)
		if not param_obj.call_method_id or not param_obj.phenotype_method_id or not param_obj.analysis_method_id:
			param_obj.results_id = request.params.get('results_id', None)
			if param_obj.results_id:
				rm = ResultsMethod.get(param_obj.results_id)
				param_obj.call_method_id = rm.call_method_id
				param_obj.phenotype_method_id = rm.phenotype_method_id
				param_obj.analysis_method_id = rm.analysis_method_id
		if param_obj.position is not None:
			param_obj.position = int(param_obj.position)
		if param_obj.call_method_id is not None:
			param_obj.call_method_id = int(param_obj.call_method_id)
		if param_obj.phenotype_method_id is not None:
			param_obj.phenotype_method_id = int(param_obj.phenotype_method_id)
		if param_obj.analysis_method_id is not None:
			param_obj.analysis_method_id = int(param_obj.analysis_method_id)
		
	def getSNPSummaryInfo(self, id=None):
		"""
		2009-6-5
			cast c.chromosome, c.position into integers. actually no need for that, pylons/elixir would cast automatically.
		2009-3-2
			information for the SNP summary table
		"""
		self.readInRequestParams(request, c)
		#3rd finally construct the full data and turn it into json
		column_name_type_ls = [("chromosome", ("string", "Chromosome")), ("position", ("string", "Position")), \
							("call_method_id", ("string", "Call Method ID")), ("phenotype", ("string","Phenotype")), \
							("analysis",("string", "Association Method")), ("score", ("number","Score")), ("type_of_snp",("string", "Type of SNP")), \
							("gene_on_the_left", ("number","Gene on the Left")), ("gene_touched", ("string", "Gene touched")), ("gene_on_the_right", ("number", "Gene on the Right"))]
		
		description = dict(column_name_type_ls)
		
		c.gene_desc_names = ['gene_symbol', 'description', 'type_of_gene', 'dbxrefs']
		
		snp = Snps.query.filter_by(chromosome=int(c.chromosome)).filter_by(position=int(c.position)).filter(Snps.end_position==None).first()
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
				gene_touched = '<a href='+h.NCBIGeneDBURL%gene.gene_id+'  target="_blank">%s %s</a>'%\
					(gene.gene_id, gene.gene_symbol)
			if snp_annotation.comment:
				snp_annotation_text += ':%s'%snp_annotation.comment
			if len(snp_annotation_text_ls)>0:
				if snp_annotation_text!=snp_annotation_text_ls[-1]:
					snp_annotation_text_ls.append(snp_annotation_text)
			else:
				snp_annotation_text_ls.append(snp_annotation_text)
		type_of_snp = '; '.join(snp_annotation_text_ls)
		
		phenotype_method = PhenotypeMethod.get(c.phenotype_method_id)
		analysis_method = AnalysisMethod.get(c.analysis_method_id)
		phenotype = '%s %s'%(phenotype_method.id, phenotype_method.short_name)
		
		analysis = '%s %s'%(analysis_method.id, analysis_method.short_name)
		
		return_ls = [dict(chromosome=c.chromosome, position=c.position, call_method_id=c.call_method_id, \
						phenotype=phenotype,\
						analysis=analysis, score=float(c.score), type_of_snp=type_of_snp, \
						gene_on_the_left=None, gene_touched=gene_touched, gene_on_the_right=None)]
		
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		return json_result
	
	def getSignificantHits(self, id=None):
		"""
		2009-6-5
			cast c.chromosome, c.position into integers. actually no need for that, pylons/elixir would cast automatically. 
		2009-4-7
			enlarge the GBrowse region from 20kb to 80kb
		2009-3-2
			significant hits (stored in Stock_250kDB.Results, assume top 1000) in all phenotypes for the given SNP
			and feed back to the javascript request
		"""
		self.readInRequestParams(request, c)
		
		column_name_type_ls = [("results_id", ("number", "Result Id")), ("phenotype_method_id", ("number", "Phenotype ID")), ("phenotype_short_name", ("string", "Phenotype Name")), \
							("analysis", ("string", "Association Method")), \
							("gbrowseLink", ("string","GBrowse")), ("score", ("number","Score")), \
							("rank", ("number","Rank")), ('beta', ('number', 'Beta')),\
							("maf", ("number","MAF")), ("mac", ("number","MAC")),\
							('genotype_var_perc', ('number', 'variance explained'))]
		
		description = dict(column_name_type_ls)
		
		snp = Snps.query.filter_by(chromosome=int(c.chromosome)).filter_by(position=int(c.position)).filter(Snps.end_position==None).first()
		
		rows = Results.query.filter_by(snps_id=snp.id).filter(Results.result.has(call_method_id=c.call_method_id)).all()
		return_ls =[]
		start_pos = c.position-40000
		stop_pos = c.position+40000
		for row in rows:
			phenotype_method = row.result.phenotype_method
			analysis_method = row.result.analysis_method
			
			analysis = '%s %s'%(analysis_method.id, analysis_method.short_name)
			track_id = '%s_%s_%s'%(row.result.call_method.id, phenotype_method.id, analysis_method.id)
			SNPURL = h.url_for(controller='SNP', action=None, phenotype_method_id=phenotype_method.id, \
									call_method_id=row.result.call_method.id, analysis_method_id=analysis_method.id,\
									chromosome=c.chromosome, position=c.position, score=row.score)
			#h.GBrowseURL%(start_pos, stop_pos, c.chromosome)+track_id+"-"+track_id+"_SNP"
			gbrowseLink = "<a href=%s target='_blank'>%s</a>"%(SNPURL, row.result.id)
			
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
	
	def getEcotypeAllelePhenotype(self, id=None):
		"""
		2009-6-22
			if phenotype_method.data_type is binary or ordered_categorical, add 1 to phenotype_value
				to make accessions with zero-phenotype-value visible in the "Ecotype Allele Phenotype BarChart"
		2009-3-5
			return json data of ecotype, allele, phenotype for the javascript in templates/snp.html
		"""
		self.readInRequestParams(request, c)
		#1st get call_info_id 2 PC12 mapping
		rows = model.Stock_250kDB.CallInfo.query.filter_by(method_id=c.call_method_id).all()
		call_info_id2pcs = {}
		for row in rows:
			call_info_id2pcs[row.id] = row.pc_value_ls[:10]
		
		#if h.ecotype_info is None:	#2009-3-6 not used right now
		#	h.ecotype_info = getEcotypeInfo(model.db)
		
		pheno_data = getattr(model, 'pheno_data', None)
		if pheno_data is None:
			from variation.src.OutputPhenotype import OutputPhenotype
			pheno_data = OutputPhenotype.getPhenotypeData(model.db.metadata.bind, phenotype_avg_table=model.Stock_250kDB.PhenotypeAvg.table.name,\
														phenotype_method_table=model.Stock_250kDB.PhenotypeMethod.table.name)
			model.pheno_data = pheno_data
		
		snpData = hc.getSNPDataGivenCallMethodID(c.call_method_id)
		pm = PhenotypeMethod.get(c.phenotype_method_id)
		
		column_name_type_ls = [("label", ("string","ID Name Phenotype")), ("date", ("date", "Date")), \
							("lon",("number", "Longitude")), ("lat",("number", "Latitude")), \
							("ecotypeid", ("number", "Ecotype ID")), ("name", ("string", "Native Name")), \
							("phenotype", ("number", "Phenotype")), \
							("allele", ("string", "Allele")), ("country", ("string", "Country")),\
							("pc1", ("number","PC1")), ("pc2", ("number", "PC2")), ("pc3", ("number","PC3")), \
							("pc4", ("number", "PC4")), ("pc5", ("number","PC5")), ("pc6", ("number", "PC6"))]
		
		description = dict(column_name_type_ls)
		rows = model.db.metadata.bind.execute("select * from view_call where call_method_id=%s"%(c.call_method_id))
		return_ls = []
		for row in rows:
			pcs = call_info_id2pcs.get(row.call_info_id)
			if pcs:
				pc_value_ls = [getattr(call_info_pc, 'pc_value') for call_info_pc in pcs]
			else:
				pc_value_ls = [0]*10
			pheno_row_index = pheno_data.row_id2row_index.get(row.ecotype_id)
			pheno_col_index = pheno_data.col_id2col_index.get(c.phenotype_method_id)
			if pheno_row_index is not None and pheno_col_index is not None:
				phenotype_value = pheno_data.data_matrix[pheno_row_index][pheno_col_index]
				if phenotype_value == 'NA':
					phenotype_value = None
			else:
				phenotype_value = None	#numpy.nan can't be recognized by ToJSon()
			label = '%s ID:%s Phenotype:%s.'%(row.nativename, row.ecotype_id, phenotype_value)
			
			snpdata_row_index = snpData.row_id2row_index.get(str(row.ecotype_id))
			snpdata_col_index = snpData.col_id2col_index.get('%s_%s'%(c.chromosome, c.position))
			if snpdata_row_index is None or snpdata_col_index is None:
				allele = -2
			else:
				allele = snpData.data_matrix[snpdata_row_index][snpdata_col_index]
			allele = number2nt[allele]
			if (pm.data_type=='binary' or pm.data_type=='ordered_categorical') and phenotype_value is not None:
				phenotype_value += 1
			return_ls.append(dict(date=datetime.date(2009,2,3), ecotypeid=row.ecotype_id, label=label, name=row.nativename, \
								lat=row.latitude, lon=row.longitude,\
								pc1=pc_value_ls[0], pc2=pc_value_ls[1], pc3=pc_value_ls[2], pc4=pc_value_ls[3], \
								pc5=pc_value_ls[4], pc6=pc_value_ls[5], \
								phenotype=phenotype_value, allele=allele, country=row.country))
		
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		#column_ls.sort()	#['date', 'ecotypeid', 'label', 'lat', 'lon', 'name', 'pc1', 'pc2', 'phenotype']
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		return json_result