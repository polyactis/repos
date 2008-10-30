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
		SNPRegionPlotType = model.Stock_250kDB.SNPRegionPlotType
		rows = SNPRegionPlotType.query.order_by(SNPRegionPlotType.id).all()
		c.snp_region_plot_type_ls = rows
		return render('/snp_region_plot.html')
	
	def getSNPInfo(self, plot_type_id=None):
		"""
		2008-10-17
			return data structure contain sorted snp list and its info
		"""
		if plot_type_id is None:
			plot_type_id = request.params.get('plot_type_id', 1)
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		rows = model.db.metadata.bind.execute("select distinct id, chromosome, start, stop, center_snp_position from %s \
			where plot_type_id=%s order by chromosome, start, center_snp_position, stop "\
			%(SNPRegionPlot.table.name, plot_type_id))
		snp_ls = []
		snp2index = {}
		snp_label_ls = []
		for row in rows:
			snp_key = (row.chromosome, row.start, row.stop)
			if snp_key not in snp2index:
				snp2index[snp_key] = len(snp_ls)
				row.no_of_phenotypes = 1
				snp_ls.append(row)
				snp_label_ls.append('%s_%s_%s'%(row.chromosome, row.start, row.stop))
			else:
				snp_index = snp2index[snp_key]
				snp_ls[snp_index].no_of_phenotypes += 1
		snp_info = PassingData()
		snp_info.snp2index = snp2index
		snp_info.snp_ls = snp_ls
		snp_info.snp_label_ls = snp_label_ls
		return snp_info
	
	def type_by_snp(self, id=None):
		"""
		2008-10-17
			display SNP region plots in a matrix fashion by spanning the phenotype axis.
		"""
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		if id is None:
			id = request.params.get('id', 1)
		extra_condition = 's.plot_type_id=%s'%id
		
		c.phenotype_info  = h.getPhenotypeInfo(SNPRegionPlot.table.name, extra_condition)
		c.snp_info = self.getSNPInfo(plot_type_id=id)
		
		data_matrix = numpy.zeros([len(c.snp_info.snp_ls), len(c.phenotype_info.phenotype_method_id_ls)], numpy.int)
		data_matrix[:] = -1
		counter = 0
		
		rows = SNPRegionPlot.query.filter_by(plot_type_id=id).all()
		for row in rows:
			snp_key = (row.chromosome, row.start, row.stop)
			row_index = c.snp_info.snp2index[snp_key]
			col_index = c.phenotype_info.phenotype_method_id2index[row.phenotype_method_id]
			data_matrix[row_index, col_index] = row.id
			counter +=1
		c.counter = counter
		c.data_matrix = data_matrix
		return render('/snp_region_plot_type_by_snp.html')
	
	def getGeneInfo(self, plot_type_id=None):
		"""
		2008-10-17
			return data structure containing sorted gene list and its info
		"""
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		SNPRegionPlotToGene = model.Stock_250kDB.SNPRegionPlotToGene
		rows = model.db.metadata.bind.execute("select g.gene_id, count(distinct s.phenotype_method_id) as count from %s g, %s s\
			where s.id=g.plot_id and s.plot_type_id=%s group by g.gene_id order by count desc, gene_id"\
			%(SNPRegionPlotToGene.table.name, SNPRegionPlot.table.name, plot_type_id))
		gene_ls = []
		gene_id2index = {}
		gene_label_ls = []
		gene_desc_names = ['gene_id', 'gene_symbol', 'type_of_gene', 'chr', 'start', 'stop', 'protein_label', 'protein_comment', 'protein_text']
		for row in rows:
			gene_id2index[row.gene_id] = len(gene_ls)
			matrix_of_gene_descriptions = h.returnGeneDescLs(model.gene_annotation, gene_id_ls=[row.gene_id])
			if len(matrix_of_gene_descriptions)>0:
				gene_desc_ls = matrix_of_gene_descriptions[0]
				gene = PassingData(gene_id=row.gene_id, gene_symbol=gene_desc_ls[1], type_of_gene=gene_desc_ls[2], chr=gene_desc_ls[3],\
								start=gene_desc_ls[4], stop=gene_desc_ls[5], protein_label=gene_desc_ls[6], protein_comment=gene_desc_ls[7], \
								protein_text=gene_desc_ls[8])
			else:
				gene = PassingData(gene_id=row.gene_id)
			gene.count = row.count
			gene.snp_region_plot_ls = []
			gene_ls.append(gene)
			gene_label_ls.append(getattr(gene, 'gene_symbol', row.gene_id))
		gene_info = PassingData()
		gene_info.gene_id2index = gene_id2index
		gene_info.gene_ls = gene_ls
		gene_info.gene_label_ls = gene_label_ls
		return gene_info
	
	def type_by_gene(self, id=None):
		"""
		2008-10-17
		"""
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		SNPRegionPlotToGene = model.Stock_250kDB.SNPRegionPlotToGene
		if id is None:
			id = request.params.get('id', 1)
		extra_condition = 's.plot_type_id=%s'%id
		
		c.phenotype_info  = h.getPhenotypeInfo(SNPRegionPlot.table.name, extra_condition)
		c.gene_info = self.getGeneInfo(plot_type_id=id)
		
		data_matrix = numpy.zeros([len(c.gene_info.gene_ls), len(c.phenotype_info.phenotype_method_id_ls)], numpy.int)
		#data_matrix[:] = -1
		counter = 0
		
		rows = SNPRegionPlotToGene.query.filter(SNPRegionPlotToGene.snp_region_plot.has(plot_type_id=id)).all()
		for row in rows:
			row_index = c.gene_info.gene_id2index[row.gene_id]
			c.gene_info.gene_ls[row_index].snp_region_plot_ls.append(row.snp_region_plot)
			col_index = c.phenotype_info.phenotype_method_id2index[row.snp_region_plot.phenotype_method_id]
			data_matrix[row_index, col_index] = row.snp_region_plot.id
			counter +=1
		for row_index in range(len(c.gene_info.gene_ls)):
			c.gene_info.gene_ls[row_index].snp_region_plot_ls.sort(lambda x, y: int(x.phenotype_method_id-y.phenotype_method_id))	#comparison function must return int. and long is NOT int.
			c.gene_info.gene_ls[row_index].snp_region_plot_id_ls = [str(snp_region_plot.id) for snp_region_plot in c.gene_info.gene_ls[row_index].snp_region_plot_ls]
		c.counter = counter
		c.data_matrix = data_matrix
		return render('/snp_region_plot_type_by_gene.html')
	
	def type(self, id=None):
		"""
		2008-10-06
		"""
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		rows = SNPRegionPlot.query.filter_by(plot_type_id=id).order_by(SNPRegionPlot.chromosome).order_by(SNPRegionPlot.start).order_by(SNPRegionPlot.phenotype_method_id).all()
		c.snp_region_plots = rows
		return render('/snp_region_plot_type.html')
	
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
				img_data = snp_region_plot.png_data
				#img_data = base64.b64decode(snp_region_plot.img_data)
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
	
	def show_plot(self, id=None, snp_region_plot_ids=None):
		"""
		2008-10-30
			add option snp_region_plot_ids for display_results_gene_one.html to call but not necessary. request.params is good enough.
		2008-10-17
			get request params: get_other_phenotypes, snp_region_plot_ids
		2008-10-06
		"""
		SNPRegionPlot = model.Stock_250kDB.SNPRegionPlot
		if id:
			c.snp_region_plot = SNPRegionPlot.get(id)
		else:
			c.snp_region_plot = None
		
		get_other_phenotypes = request.params.get('get_other_phenotypes', False)
		snp_region_plot_ids = request.params.get('snp_region_plot_ids', snp_region_plot_ids)
		if get_other_phenotypes:
			c.snp_region_plot_ls = SNPRegionPlot.query.filter_by(chromosome=c.snp_region_plot.chromosome).\
				filter_by(start=c.snp_region_plot.start).\
				filter_by(stop=c.snp_region_plot.stop).\
				filter_by(plot_type_id=c.snp_region_plot.plot_type_id).order_by(SNPRegionPlot.phenotype_method_id).all()
		elif snp_region_plot_ids:
			snp_region_plot_id_ls = snp_region_plot_ids.split(',')
			#snp_region_plot_id_ls = map(int, snp_region_plot_id_ls)
			c.snp_region_plot_ls = [SNPRegionPlot.get(int(snp_region_plot_id)) for snp_region_plot_id in snp_region_plot_id_ls]
		else:
			c.snp_region_plot_ls = []
		
		if c.snp_region_plot is None and len(c.snp_region_plot_ls)>0:
			c.snp_region_plot = c.snp_region_plot_ls[0]
		
		"""
		if 'gene_annotation' not in session:
			session['gene_annotation'] = h.dealWithGeneAnnotation()
			session.save()
		"""
		#plot_file_path = os.path.join(config['app_conf']['plots_store'], 'snp_region_plot_%s.png'%c.snp_region_plot.id)
		#if not os.path.isfile(plot_file_path):
		#	self.savePlot(c.snp_region_plot, plot_file_path)
		
		gene_annotation = model.gene_annotation
		gene_id_ls = [plot2gene.gene_id for plot2gene in c.snp_region_plot.plot2gene_ls]
		c.matrix_of_gene_descriptions = h.returnGeneDescLs(gene_annotation, gene_id_ls)
		#print abc
		#raise Exception('Just testing the interactive debugger!')
		#c.matrix_of_gene_descriptions = []
		return render('/snp_region_plot_single.html')