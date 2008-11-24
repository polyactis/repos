#!/usr/bin/env python
"""
Examples:
	DrawFTPathway.py -o /tmp/phenotype_1_a1_m5000_g0 -y 1 -a 1 -m 5000
	
	#draw flowering pathways for all FT phenotypes with method 1(KW) and 7(emma)
	for y in 1 2 3 4 5 6 7 39 40 41 42 43 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 80 81 82; do echo $y; for a in 1 7 ; do echo $a; ./script/variation/src/DrawFTPathway.py -o /Network/Data/250k/tmp-yh/FTPathway/phenotype_$y\_a$a\_m5000_g0 -y $y -a $a  -u yh -p passw***; done ; done
	
Description:
	Program to draw Flowering time pathway with genes colored according to pvalue significance. Data is all in db, including ft pathway genes.
"""

import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib as mpl; mpl.use("Agg")	#to avoid popup and collapse in X11-disabled environment
from matplotlib import rcParams
from PlotGroupOfSNPs import PlotGroupOfSNPs
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 6
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 6
rcParams['axes.titlesize'] = 9
rcParams['xtick.labelsize'] = 4
rcParams['ytick.labelsize'] = 4

import time, csv, cPickle
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr, SNPData, read_data,\
	assignMatPlotlibHueColorToLs, drawName2FCLegend
from pymodule.DrawMatrix import Value2Color
import Stock_250kDB, StockDB
from sets import Set
from GeneListRankTest import GeneListRankTest	#GeneListRankTest.getGeneList()

from matplotlib.patches import Polygon, CirclePolygon, Ellipse, Wedge
import pylab
import ImageColor
import numpy
from Kruskal_Wallis import Kruskal_Wallis
from PhenotypeOfAncestralDerivedAllele import PhenotypeOfAncestralDerivedAllele
from common import get_chr_id2size, get_chr_id2cumu_size, getEcotypeInfo

class DrawFTPathway(PlotGroupOfSNPs):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('output_fname_prefix', 1, ): ['', 'o', 1, 'store the pvalue', ],\
							('phenotype_method_id', 1, int): [None, 'y', 1, 'which phenotype',],\
							('min_MAF', 1, float): [0, 'n', 1, 'minimum Minor Allele Frequency.'],\
							('min_distance', 1, int): [5000, 'm', 1, 'minimum distance'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							('call_method_id', 0, int):[17, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id', 0, int):[7, 'a', 1, 'Restrict results based on this analysis_method. Default is no such restriction.'],\
							('score_is_rank', 0, ): [0, '', 0, 'take rank as score, rather than -log(pvalue)', ],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'debug mode. 1=level 1 (pdb mode). 2=level 2 (same as 1 except no pdb mode)'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-11-14
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def get_results_by_gene(self, call_method_id, analysis_method_id, phenotype_method_id, \
			min_distance, get_closest, min_MAF, results_directory, score_is_rank=0):
		"""
		"""
		sys.stderr.write("Getting results_by_gene ...")
		rm = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).\
					filter_by(phenotype_method_id=phenotype_method_id).filter_by(analysis_method_id=analysis_method_id).first()
		ResultsByGene = Stock_250kDB.ResultsByGene
		rows = ResultsByGene.query.filter_by(results_method_id=rm.id).filter_by(min_distance=min_distance).\
			filter_by(get_closest=get_closest).\
			filter(ResultsByGene.min_MAF>=min_MAF-0.0001).filter(ResultsByGene.min_MAF<=min_MAF+0.0001)
		rbg = rows.first()
		if not rbg:
			sys.stderr.write("No results_by_gene available.\n")
			return None
		
		
		if results_directory:	#given a directory where all results are.
			result_fname = os.path.join(results_directory, os.path.basename(rbg.filename))
		else:
			result_fname = rbg.filename
		reader = csv.reader(open(result_fname), delimiter='\t')
		reader.next()
		gene_id2score = {}
		counter = 0
		for row in reader:
			counter += 1
			gene_id = int(row[0])
			if score_is_rank:
				score = -float(counter)	#lower rank is more significant. in order to keep the correlation between score and significance, minus rank.
			else:
				score = float(row[1])
			gene_id2score[gene_id] = score
		sys.stderr.write("Done.\n")
		return gene_id2score
	
	def get_pathway_id2gene_id_ls(self):
		"""
		2008-11-14
		"""
		sys.stderr.write("Getting pathway_id2gene_id_ls ... ")
		pathway_id2gene_id_ls = {}
		rows = Stock_250kDB.FTGene.query.all()
		for row in rows:
			if row.pathway_id not in pathway_id2gene_id_ls:
				pathway_id2gene_id_ls[row.pathway_id] = []
			pathway_id2gene_id_ls[row.pathway_id].append(row.gene_id)
			
		sys.stderr.write("Done.\n")
		return pathway_id2gene_id_ls
	
	def get_min_max_score(self, pathway_id2gene_id_ls, gene_id2score):
		"""
		2008-11-15
		"""
		min_score =None
		max_score = None
		for pathway_id, gene_id_ls in pathway_id2gene_id_ls.iteritems():
			for gene_id in gene_id_ls:
				score = gene_id2score.get(gene_id)
				if score:
					if min_score==None:
						min_score = score
					elif score<min_score:
						min_score = score
					
					if max_score == None:
						max_score = score
					elif score>max_score:
						max_score = score
		return min_score , max_score
	
	def return_rgb_color_given_score(self, score, value2color_func):
		"""
		2008-11-15
		"""
		hue_value_in_hsl_format = value2color_func(score)
		fc = ImageColor.getrgb(hue_value_in_hsl_format)	#"hsl(hue, saturation%, lightness%)" where hue is the colour given as an
		# angle between 0 and 360 (red=0, green=120, blue=240),
		#saturation is a value between 0% and 100% (gray=0%, full color=100%), and lightness is a value between 0% and 100% (black=0%, normal=50%, white=100%).
		fc = [color_value/255. for color_value in fc]
		return fc
	
	def draw_pathway_box(self, gene_id_ls, axe, gene_id2gene_symbol, gene_id2score, value2color_func, x_span=1., y_span=1., \
			gene_label_size=8):
		"""
		"""
		sys.stderr.write("drawing pathway box ...")
		no_of_genes_per_line = min(4, len(gene_id_ls)) 
		gene_label_width = x_span/no_of_genes_per_line
		no_of_lines = math.ceil(len(gene_id_ls)/float(no_of_genes_per_line))
		no_of_lines = int(no_of_lines)
		gene_label_height = y_span/no_of_lines
		
		for i in range(len(gene_id_ls)):
			gene_id = gene_id_ls[i]
			which_line = i/no_of_genes_per_line
			which_col = i%no_of_genes_per_line
			x_pos = gene_label_width/2+which_col*gene_label_width
			y_pos = 1-(gene_label_height/2 + which_line*gene_label_height)
			score = gene_id2score[gene_id]
			gene_symbol = gene_id2gene_symbol[gene_id]
			
			fc = self.return_rgb_color_given_score(score, value2color_func)
			center = (x_pos, y_pos)
			patch = Ellipse(center, gene_label_width*0.75, gene_label_height/2, facecolor=fc, linewidth=0)
			axe.add_patch(patch)
			
			axe.text(x_pos, y_pos, gene_symbol, \
								horizontalalignment ='center', verticalalignment='center', size=gene_label_size)
			
		sys.stderr.write("Done.\n")
	
	def return_xy_given_loc(self, arrow_start_point_loc='bottom middle'):
		arrow_start_point_loc = arrow_start_point_loc.split(' ')
		if arrow_start_point_loc[0]=='bottom':
			start_y_pos = 0
			if arrow_start_point_loc[1] =='left':
				start_x_pos = 0
			elif arrow_start_point_loc[1] =='right':
				start_x_pos = 1
			elif arrow_start_point_loc[1] == 'middle':
				start_x_pos = 0.5
			else:
				start_x_pos = 0.5	
		elif arrow_start_point_loc[0] == 'top':
			start_y_pos = 1
			if arrow_start_point_loc[1] =='left':
				start_x_pos = 0
			elif arrow_start_point_loc[1] =='right':
				start_x_pos = 1
			elif arrow_start_point_loc[1] == 'middle':
				start_x_pos = 0.5
			else:
				start_x_pos = 0.5
		elif arrow_start_point_loc[0]=='left':
			start_x_pos = 0.
			if arrow_start_point_loc[1] =='top':
				start_y_pos = 1
			elif arrow_start_point_loc[1] =='bottom':
				start_y_pos = 0
			elif arrow_start_point_loc[1] == 'middle':
				start_y_pos = 0.5
			else:
				start_y_pos = 0.5
		elif arrow_start_point_loc[0]=='right':
			start_x_pos = 1.
			if arrow_start_point_loc[1] =='top':
				start_y_pos = 1
			elif arrow_start_point_loc[1] =='bottom':
				start_y_pos = 0
			elif arrow_start_point_loc[1] == 'middle':
				start_y_pos = 0.5
			else:
				start_y_pos = 0.5
		return (start_x_pos, start_y_pos)
	
	def draw_edge_between_pathways(self, axe_cover, axe_start, axe_stop, edge_type=1, arrow_start_point_loc ='bottom middle',\
			arrow_end_point_loc='top middle'):
		"""
		2008-11-15
		"""
		sys.stderr.write("drawing edges between pathways ...")
		start_x_pos, start_y_pos = self.return_xy_given_loc(arrow_start_point_loc)
		stop_x_pos, stop_y_pos = self.return_xy_given_loc(arrow_end_point_loc)
		
		#freeze transformations of all 3 axes
		for ax in [axe_cover, axe_start, axe_stop]:
			ax.transData.freeze()  # eval the lazy objects
			ax.transAxes.freeze()
		start_x_canvas, start_y_canvas = axe_start.transData.xy_tup((start_x_pos, start_y_pos))
		start_x_cover,  start_y_cover = axe_cover.transData.inverse_xy_tup( (start_x_canvas, start_y_canvas) )
		
		stop_x_canvas, stop_y_canvas = axe_stop.transData.xy_tup((stop_x_pos, stop_y_pos))
		stop_x_cover,  stop_y_cover = axe_cover.transData.inverse_xy_tup( (stop_x_canvas, stop_y_canvas))
		if edge_type == 1:
			arrow_color = 'g'
		elif edge_type==-1:
			arrow_color = 'r'
		else:
			arrow_color = 'k'
		axe_cover.plot([start_x_cover, stop_x_cover] , [start_y_cover, stop_y_cover], c=arrow_color, \
												linestyle='-', alpha=0.8, linewidth=2, zorder=0)
		
		for ax in [axe_cover, axe_start, axe_stop]:
			ax.transData.thaw() 
			ax.transAxes.thaw()
		sys.stderr.write("Done.\n")
		
	def drawLegend(self, ax, value2color_func, min_score, max_score, no_of_bands =20, label_size=8, score_is_rank=True):
		"""
		2008-11-14
		"""
		sys.stderr.write("\t Drawing legend  ...")
		xs = [0,0,0.5,0.5]	#X-axis value for the 4 points of the rectangle starting from upper left corner. 
		ys = numpy.array([0.1,0,0,0.1])
		score_gap  = (max_score-min_score)/float(no_of_bands)
		ys_gap = 1.0/no_of_bands
		for i in range(no_of_bands):
			score = min_score + i*score_gap
			fc = self.return_rgb_color_given_score(score, value2color_func)
			poly = Polygon(zip(xs, ys), facecolor=fc, linewidth=0)
			ax.add_patch(poly)
			if i%3==0:	#every other band
				if score_is_rank:
					score_label = int(abs(score))
				else:
					score_label = '%.2f'%score
				ax.text(0.6, ys[0], score_label, horizontalalignment ='left', verticalalignment='top', size=label_size)
			ys += ys_gap	#increase y-axis
		if score_is_rank:
			max_score_label = "rank=%s"%int(abs(max_score))
		else:
			max_score_label = 'score=%.2f'%max_score
		ax.text(0.6, 1.1, max_score_label, horizontalalignment ='left', verticalalignment='top', size=label_size)
		sys.stderr.write("Done.\n")
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		import MySQLdb
		mysql_conn = MySQLdb.connect(db=self.dbname, host='banyan.usc.edu', user = self.db_user, passwd = self.db_passwd)
		mysql_curs = mysql_conn.cursor()
		from pymodule.utils import get_gene_id2gene_symbol
		gene_id2gene_symbol = get_gene_id2gene_symbol(mysql_curs, 3702, table='genome.gene', upper_case_gene_symbol=0)	#3702 is At's tax id
		
		self.get_closest = 0
		gene_id2score = self.get_results_by_gene(self.call_method_id, self.analysis_method_id, self.phenotype_method_id,\
				self.min_distance, self.get_closest, self.min_MAF, self.results_directory, self.score_is_rank)
		
		pathway_id2gene_id_ls = self.get_pathway_id2gene_id_ls()
		
		min_score, max_score = self.get_min_max_score(pathway_id2gene_id_ls, gene_id2score)
		
		value2color_func = lambda x: Value2Color.value2HSLcolor(x, min_score, max_score)
		
		
		if not gene_id2score:
			sys.stderr.write("gene_id2score is empty. exit.\n")
			sys.exit(3)
		
		gap = 0.05
		
		axe_auto_x_offset = 0.05
		axe_auto_y_offset = 0.7
		axe_auto_width = 0.3
		axe_auto_height = 0.2
		axe_auto = pylab.axes([axe_auto_x_offset, axe_auto_y_offset, axe_auto_width, axe_auto_height])
		
		
		axe_vern_x_offset = axe_auto_x_offset + axe_auto_width + gap
		axe_vern_y_offset = axe_auto_y_offset
		axe_vern_width = axe_auto_width*0.75
		axe_vern_height = axe_auto_height
		axe_vern = pylab.axes([axe_vern_x_offset, axe_vern_y_offset, axe_vern_width, axe_vern_height])
		
		axe_light_x_offset = axe_vern_x_offset + axe_vern_width + gap
		axe_light_y_offset = axe_vern_y_offset
		axe_light_width = axe_auto_width
		axe_light_height = axe_auto_height
		axe_light = pylab.axes([axe_light_x_offset, axe_light_y_offset, axe_light_width, axe_light_height])
		
		axe_flc_act_x_offset = 0.05
		axe_flc_act_y_offset = 0.45
		axe_flc_act_width = 0.3
		axe_flc_act_height = 0.2
		axe_flc_act = pylab.axes([axe_flc_act_x_offset, axe_flc_act_y_offset, axe_flc_act_width, axe_flc_act_height])

		axe_repr_x_offset = 0.05
		axe_repr_y_offset = 0.1
		axe_repr_width  = 0.3
		axe_repr_height = 0.2
		axe_repr = pylab.axes([axe_repr_x_offset, axe_repr_y_offset, axe_repr_width, axe_repr_height] )
		
		axe_flc_x_offset = 0.4
		axe_flc_y_offset = 0.5
		axe_flc_width = 0.1
		axe_flc_height = axe_flc_act_height/3
		axe_flc = pylab.axes([axe_flc_x_offset, axe_flc_y_offset, axe_flc_width, axe_flc_height])
		
		axe_co_x_offset = 0.6
		axe_co_y_offset = axe_flc_y_offset
		axe_co_width = axe_flc_width
		axe_co_height = axe_flc_height
		axe_co = pylab.axes([axe_co_x_offset, axe_co_y_offset, axe_co_width, axe_co_height])
		
		axe_gibbe_x_offset = axe_co_x_offset + axe_co_width + gap
		axe_gibbe_y_offset = axe_co_y_offset
		axe_gibbe_width = axe_auto_width*0.75
		axe_gibbe_height = axe_co_height
		axe_gibbe = pylab.axes([axe_gibbe_x_offset, axe_gibbe_y_offset, axe_gibbe_width, axe_gibbe_height ])
		
		
		axe_int_x_offset = 0.4
		axe_int_y_offset = 0.3
		axe_int_width = 0.3
		axe_int_height = 0.05
		axe_int = pylab.axes([axe_int_x_offset, axe_int_y_offset, axe_int_width, axe_int_height])
		
		axe_meri_x_offset = 0.4
		axe_meri_y_offset = 0.1
		axe_meri_width = axe_int_width
		axe_meri_height = axe_int_height
		axe_meri = pylab.axes([axe_meri_x_offset, axe_meri_y_offset, axe_meri_width, axe_meri_height])
		
		axe_flower = pylab.axes([axe_meri_x_offset, 0.01, axe_meri_width, 0.05])
		axe_flower.text(0.5, 0.5, 'Flowering', \
								horizontalalignment ='center', verticalalignment='center', size=8)
		#iteration over pathway_id2axe will do. omit this section.
		axe_flower.set_xticks([])
		axe_flower.set_yticks([])
		
		axe_score_legend = pylab.axes([axe_int_x_offset+axe_int_width+0.1, 0.01, 0.1, 0.3], frameon=False)	#axes for the legend of LD
		axe_score_legend.set_xticks([])
		axe_score_legend.set_yticks([])
		
		axe_cover = pylab.axes([0.03,0.01, 0.95,0.95],frameon=False)
		axe_cover.set_xticks([])
		axe_cover.set_yticks([])
		
		pathway_id2axe = {1:axe_flc_act, 2: axe_repr, 3:axe_light, 4:axe_gibbe, 5:axe_vern, 6:axe_int, 7: axe_meri, \
						8:axe_flc, 9:axe_co, 10:axe_auto}
		#draw edges first to avoid the edges overwriting the pathway labels
		pathway_edges = Stock_250kDB.FTPathwayRelationship.query.all()
		
		for pathway_edge in pathway_edges:
			
			axe_start = pathway_id2axe[pathway_edge.pathway1_id]
			axe_stop = pathway_id2axe[pathway_edge.pathway2_id]
			axe_cover.set_xlim([0,1])
			axe_cover.set_ylim([0,1])
			self.draw_edge_between_pathways(axe_cover, axe_start, axe_stop, edge_type=pathway_edge.relationship_type_id, \
				arrow_start_point_loc =pathway_edge.arrow_start_point_loc,\
				arrow_end_point_loc=pathway_edge.arrow_end_point_loc)
			axe_cover.set_xlim([0,1])
			axe_cover.set_ylim([0,1])
		
		for pathway_id, axe in pathway_id2axe.iteritems():
			axe.set_xticks([])
			axe.set_yticks([])
			pathway = Stock_250kDB.FTPathway.get(pathway_id)
			axe.title.set_text(pathway.short_name)
			gene_id_ls = pathway_id2gene_id_ls.get(pathway_id)
			if gene_id_ls:
				self.draw_pathway_box(gene_id_ls, axe, gene_id2gene_symbol, gene_id2score, value2color_func, x_span=1., y_span=1., \
						gene_label_size=8)
		
		#draw the final edge from axe_meri to axe_flower
		self.draw_edge_between_pathways(axe_cover, axe_meri, axe_flower, edge_type=0, arrow_start_point_loc='bottom middle',\
					arrow_end_point_loc='top middle')
		
		phenotype_method = Stock_250kDB.PhenotypeMethod.get(self.phenotype_method_id)
		analysis_method = Stock_250kDB.AnalysisMethod.get(self.analysis_method_id)
		axe_cover.title.set_text('%s by %s'%(phenotype_method.short_name, analysis_method.short_name))
		axe_cover.set_xlim([0,1])
		axe_cover.set_ylim([0,1])

		self.drawLegend(axe_score_legend, value2color_func, min_score, max_score, no_of_bands = 20, score_is_rank=self.score_is_rank)
		
		
		
		png_output_fname = '%s.png'%self.output_fname_prefix
		pylab.savefig(png_output_fname, dpi=600)
		pylab.savefig('%s.svg'%self.output_fname_prefix)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DrawFTPathway
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
