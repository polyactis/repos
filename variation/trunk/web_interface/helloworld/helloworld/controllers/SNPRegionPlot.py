import logging

from helloworld.lib.base import *
import base64, os, sys
log = logging.getLogger(__name__)

sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import helloworld.model as model
from variation.src.DrawSNPRegion import DrawSNPRegion

class SnpregionplotController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		rows = SNPRegionPlot.query.order_by(SNPRegionPlot.chromosome).order_by(SNPRegionPlot.start).order_by(SNPRegionPlot.phenotype_method_id).all()
		c.snp_region_plots = rows
		return render('/snp_region_plot.html')
	
	def type(self, id=None):
		"""
		2008-10-06
		"""
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		rows = SNPRegionPlot.query.filter_by(plot_type_id=id).order_by(SNPRegionPlot.chromosome).order_by(SNPRegionPlot.start).order_by(SNPRegionPlot.phenotype_method_id).all()
		c.snp_region_plots = rows
		return render('/snp_region_plot.html')
	
	def getImage(self, id=None):
		"""
		2008-10-06
		"""
		response.headers['Content-type'] = 'image/png'
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		if id is None:
			id = request.params.get('id', None)
		get_img_data_success = 0
		if id:
			snp_region_plot = SNPRegionPlot.get(id)
			if snp_region_plot:
				img_data = base64.b64decode(snp_region_plot.img_data)
				get_img_data_success = 1
		if not get_img_data_success:
			response.headers['Content-type'] = 'text/plain'
			img_data = 'No Image'
		return img_data
	
	def savePlot(self, snp_region_plot, image_file_path=None):
		"""
		2008-10-06
		"""
		#response.headers['Content-type'] = 'image/png'
		outf = open(image_file_path, 'wb')
		outf.write(base64.b64decode(snp_region_plot.img_data))
		outf.close()
	
	def show_plot(self, id=None):
		"""
		2008-10-06
		"""
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		c.snp_region_plot = SNPRegionPlot.get(id)
		"""
		if 'gene_annotation' not in session:
			session['gene_annotation'] = h.dealWithGeneAnnotation()
			session.save()
		"""
		#plot_file_path = os.path.join(config['app_conf']['plots_store'], 'snp_region_plot_%s.png'%c.snp_region_plot.id)
		#if not os.path.isfile(plot_file_path):
		#	self.savePlot(c.snp_region_plot, plot_file_path)
		
		gene_annotation = model.gene_annotation
		c.matrix_of_gene_descriptions = h.returnGeneDescLs(gene_annotation, c.snp_region_plot.plot2gene_ls)
		#print abc
		#raise Exception('Just testing the interactive debugger!')
		#c.matrix_of_gene_descriptions = []
		return render('/snp_region_plot_single.html')