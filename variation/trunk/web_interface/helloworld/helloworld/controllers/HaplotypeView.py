import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render, h
from helloworld import model

from DisplayResults import DisplayresultsController
from DisplayResultsGene import DisplayresultsgeneController
from variation.src.common import getEcotypeInfo
from variation.src.GeneListRankTest import GeneListRankTest
from variation.src.DrawSNPRegion import DrawSNPRegion, SNPPassingData
from HelpOtherControllers import HelpothercontrollersController as hc

log = logging.getLogger(__name__)

class HaplotypeviewController(BaseController):

	def index(self):
		#c.call_method_ls = DisplayresultsController.getCallMethodLsJson()		
		c.callMethodLsURL = h.url_for(controller="DisplayResults", action="getCallMethodLsJson", id=None)
		#c.gene_list_ls = DisplayresultsgeneController.getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethodJson()
		c.geneListLsURL = h.url_for(controller="DisplayResultsGene", action="getGeneListTypeLsGivenTypeAndPhenotypeMethodAndAnalysisMethodJson", id=None)
		c.callMethodOnChangeURL = h.url_for(controller="DisplayResults", action="getPhenotypeMethodLsJson", id=None)
		c.haplotypeImgURL = h.url_for(controller='HaplotypeView', action='getPlot', id=None)
		return render('/HaplotypeView.html')
	
	def getPlot(self):
		"""
		2009-4-30
		"""
		chromosome = int(request.params.get('chromosome', 1))
		start = int(request.params.get('start', 1))
		stop = int(request.params.get('stop', 10))
		if start>stop:
			start, stop = stop, start
		snps_id = '%s_%s'%(chromosome, start)
		
		call_method_id = int(request.params.getone('call_method_id'))
		phenotype_method_id = int(request.params.getone('phenotype_method_id'))
		list_type_id = int(request.params.get('list_type_id', None))
		
		plot_key = (chromosome, start, stop, call_method_id, phenotype_method_id, list_type_id)
		
		if not hasattr(model, "plot_key2png_data"):
			model.plot_key2png_data = {}
		
		png_data = model.plot_key2png_data.get(plot_key)
		if png_data is not None:
			if hasattr(png_data, 'getvalue'):	# png_data is StringIO
				response.headers['Content-type'] = 'image/png'
				return png_data.getvalue()
			else:
				return png_data
		
		if not getattr(model, 'phenotype_id2analysis_method_id2gwr', None):
			model.phenotype_id2analysis_method_id2gwr = {}
		analysis_method_id2gwr = model.phenotype_id2analysis_method_id2gwr.get(phenotype_method_id)
		if not analysis_method_id2gwr:
			analysis_method_id2gwr = DrawSNPRegion.getSimilarGWResultsGivenResultsByGene(phenotype_method_id, call_method_id)
			model.phenotype_id2analysis_method_id2gwr[phenotype_method_id] = analysis_method_id2gwr
		
		if list_type_id>0:	#2009-2-22
			candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
		else:
			candidate_gene_set = set()
		gene_annotation = model.gene_annotation
		
		snpData = hc.getSNPDataGivenCallMethodID(call_method_id)
		pheno_data = hc.getPhenotypeDataInSNPDataOrder(snpData)
	
		if h.ecotype_info is None:	#2009-3-6 not used right now
			h.ecotype_info = getEcotypeInfo(model.db)
		
		snp_info = getattr(model, 'snp_info', None)
		if snp_info is None:
			snp_info = DrawSNPRegion.getSNPInfo(model.db)
			model.snp_info = snp_info

		LD_info = None
		output_dir = '/tmp/'	#not used at all, place holder
		which_LD_statistic = 1
		this_snp = SNPPassingData(chromosome=chromosome, position=start, stop=stop, snps_id='%s_%s'%(chromosome, start))
		
		DrawSNPRegion.construct_chr_pos2index_forSNPData(snpData)	#prerequisite

		after_plot_data = DrawSNPRegion.drawRegionAroundThisSNP(phenotype_method_id, this_snp, candidate_gene_set, \
													gene_annotation, snp_info, \
								analysis_method_id2gwr, LD_info, output_dir, which_LD_statistic, \
								min_distance=20000, list_type_id=list_type_id,
								label_gene=True, \
								draw_LD_relative_to_center_SNP=False,\
								commit=True, snpData=snpData, phenData=pheno_data, \
								ecotype_info=h.ecotype_info, snpData_before_impute=None,\
								snp_matrix_data_type=1)
		
		if getattr(after_plot_data, 'png_data', None):
			response.headers['Content-type'] = 'image/png'
			model.plot_key2png_data[plot_key] = after_plot_data.png_data
			return after_plot_data.png_data.getvalue()
		else:
			model.plot_key2png_data[plot_key] = "No Plot"
			return model.plot_key2png_data[plot_key]