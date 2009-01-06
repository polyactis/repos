#!/usr/bin/env python
"""
2008-01-30
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

if __name__ == '__main__':
	import pygtk
	pygtk.require('2.0')
	import gtk, gtk.glade, gobject
	from gtk import gdk
	import gnome
	import gnome.ui
	import gnomecanvas
	
	import matplotlib
	matplotlib.use('GTKAgg')  # or 'GTK'
	from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
	
	from matplotlib.figure import Figure
	
	#from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
	#2008-02-04 use a custom navigation tool bar
	from pymodule.yh_gnome import NavigationToolbar2GTKAgg_chromosome as NavigationToolbar
	
	from pymodule import yh_gnome
	from variation.src.common import get_chr_pos_from_x_axis_pos

import numpy, traceback
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle, Polygon
from matplotlib.text import Text
from matplotlib.collections import LineCollection, Collection

from pymodule.yh_matplotlib_artists import Gene, ExonIntronCollection
from pymodule.db import TableClass
from Results2DB_250k import Results2DB_250k
from pymodule import GenomeWideResults, GenomeWideResult, DataObject, getGenomeWideResultFromFile, PassingData
from DrawSNPRegion import DrawSNPRegion	#2008-12-16 dealWithGeneAnnotation()
import Stock_250kDB
from GeneListRankTest import GeneListRankTest

class GeneModel:
	def __init__(self, gene_id=None, chromosome=None, symbol = None, description = None, type_of_gene = None, \
				start=None, stop=None, mrna_start = None, mrna_stop = None, cds_start = None, cds_stop = None, \
				strand = None, go_id_ls=None, go_evidence=None, go_description=None, go_term_type=None):
		self.gene_id = gene_id
		self.chromosome = chromosome
		self.symbol = symbol
		self.description = description
		self.type_of_gene = type_of_gene
		self.start = start
		self.stop = stop
		self.mrna_start = mrna_start
		self.mrna_stop = mrna_stop
		self.cds_start = cds_start
		self.cds_stop = cds_stop
		self.strand = strand
		self.go_id_ls = go_id_ls
		self.go_evidence = go_evidence
		self.go_description = go_description
		self.go_term_type = go_term_type


class GenomeBrowser(object):
	def __init__(self):
		"""
		2008-01-30
		"""
		program_path = os.path.dirname(sys.argv[0])
		xml = gtk.glade.XML(os.path.join(program_path, 'GenomeBrowser.glade'))
		xml.signal_autoconnect(self)
		self.xml = xml
		
		self.app1 = xml.get_widget("app1")
		self.app1.connect("delete_event", gtk.main_quit)
		#self.app1.set_default_size(1200, 800)
		
		self.vbox_matplotlib = xml.get_widget('vbox_matplotlib')
		
		# matplotlib canvas
		fig = Figure(figsize=(8,8))
		self.canvas_matplotlib = FigureCanvas(fig)  # a gtk.DrawingArea
		self.canvas_matplotlib.set_size_request(600,400)
		self.canvas_matplotlib.mpl_connect('pick_event', self.on_canvas_pick)
		self.vbox_matplotlib.pack_start(self.canvas_matplotlib)
		
		#matplotlib axes
		#self.ax = fig.add_subplot(111)
		axe_y_offset1 = 0.05	#y_offset for axe_LD, axe_strain_pca, axe_phenotype, axe_map
		axe_height1 = 0.15	#height of axe_LD or axe_snp_matrix
		axe_y_offset2 = axe_y_offset1+axe_height1
		axe_height2 = 0.75	#height of axe_gene_model
		axe_y_offset3 = axe_y_offset2+axe_height2
		
		
		axe_x_offset1 = 0.1	#
		axe_width1 = 0.85
		axe_x_offset2 = axe_x_offset1 + axe_width1
		
		self.ax = fig.add_axes([axe_x_offset1, axe_y_offset2, axe_width1, axe_height2], frameon=False)
		self.ax.grid(True, alpha=0.3)
		#self.ax.set_xticklabels([])	#remove xtick labels on ax1 because axe_LD's xtick labels cover this.
		
		self.axe_gene_model = fig.add_axes([axe_x_offset1, axe_y_offset1, axe_width1, axe_height1], frameon=False, sharex=self.ax)
		#axe_gene_model.set_xticks([])	#this will set ax1's xticks off as well because the x-axis is shared.
		self.axe_gene_model.set_yticks([])
		
		
		# matplotlib toolbar
		self.toolbar = NavigationToolbar(self.canvas_matplotlib, self.app1)
		self.vbox_matplotlib.pack_start(self.toolbar, False, False)
		
		self.textview_output = xml.get_widget('textview_output')
		
		self.textbuffer_output = self.textview_output.get_buffer()
		
		#redirect stdout/stderr to textbuffer_output
		t_table=self.textbuffer_output.get_tag_table()
		tag_err=gtk.TextTag("error")
		tag_err.set_property("foreground","red")
		#tag_err.set_property("font","monospace 10")
		t_table.add(tag_err)
		tag_out=gtk.TextTag("output")
		tag_out.set_property("foreground","blue")
		#tag_out.set_property("font","monospace 10")
		t_table.add(tag_out)
		
		self.dummy_out = yh_gnome.Dummy_File(self.textbuffer_output, tag_out)
		self.dummy_err = yh_gnome.Dummy_File(self.textbuffer_output, tag_err)
		sys.stdout = self.dummy_out
		sys.stderr = self.dummy_err
		
		self.app1.show_all()
		
		self.filechooserdialog1 = xml.get_widget("filechooserdialog1")
		self.entry_min_value_cutoff = xml.get_widget('entry_min_value_cutoff')
		self.filechooserdialog1.connect("delete_event", yh_gnome.subwindow_hide)
		
		
		self.dialog_db_connect = xml.get_widget("dialog_db_connect")
		self.dialog_db_connect.connect("delete_event", yh_gnome.subwindow_hide)
		self.entry_mysql_hostname = xml.get_widget("entry_mysql_hostname")
		self.entry_mysql_dbname = xml.get_widget("entry_mysql_dbname")
		self.entry_postgres_hostname = xml.get_widget("entry_postgres_hostname")
		self.entry_postgres_dbname = xml.get_widget("entry_postgres_dbname")
		self.entry_postgres_schema = xml.get_widget("entry_postgres_schema")
		self.entry_gene_annotation_picklef = xml.get_widget("entry_gene_annotation_picklef")
		self.filechooserbutton_gene_annot = xml.get_widget("filechooserbutton_gene_annot")
		
		self.dialog_preferences = xml.get_widget("dialog_preferences")
		self.dialog_preferences.connect("delete_event", yh_gnome.subwindow_hide)
		self.checkbutton_debug = xml.get_widget("checkbutton_debug")
		self.checkbutton_stdout = xml.get_widget("checkbutton_stdout")
		self.checkbutton_stderr = xml.get_widget("checkbutton_stderr")
		self.entry_gene_width = xml.get_widget("entry_gene_width")
		self.checkbutton_draw_gene_symbol = xml.get_widget("checkbutton_draw_gene_symbol")
		
		self.aboutdialog1 = xml.get_widget("aboutdialog1")
		self.aboutdialog1.connect("delete_event", yh_gnome.subwindow_hide)
		
		self.mysql_conn = self.mysql_curs = self.postgres_conn = self.postgres_curs = self.db = None
		self.gene_annotation = None
		self.candidate_gene_set = None
		
		self.chr_id2size = None
		self.chr_id2cumu_size = None
		self.chr_gap = None
		self.chr_id_ls = []
		
		self.genome_wide_results = GenomeWideResults(gap=1.0)
		self.genome_wide_results.genome_wide_result_ls = []
		self.genome_wide_results.genome_wide_result_obj_id2index = {}
		self.artist_obj_id2data_obj_key = {}
		self.yticks = []
		self.yticklabels = []
		
		self.gene_id2artist_object_id = {}
		self.chr_id2gene_id_ls = {}	#chr_id here is str type (db is varchar type)
		self.gene_id2model = {}
		self.artist_obj_id2artist_gene_id_ls = {}
		
		self.gene_id2vspan_obj_id = {}	#for the axvspan()'s drawn on the canvas
		
		self.gene_width = 1.0
		
		self.draw_gene_symbol_when_clicked = 0
		
		self.debug = 0
	
	def load_data(self, mysql_curs, postgres_curs):
		"""
		2008-02-04 update the info related to chromosome , position in toolbar
		2008-02-01
			read the input data
			
			chr_id2cumu_size has an extra fake chromosome (0) compared to chr_id2size
			
			chr_id is all changed into str type
		"""		
		from variation.src.common import get_chr_id2size, get_chr_id2cumu_size
		chr_id_int2size = get_chr_id2size(mysql_curs)
		self.chr_id2size = {}	#change the type of chr_id into string type
		for chr_id_int, size in chr_id_int2size.iteritems():
			self.chr_id2size[str(chr_id_int)] = size
		#chr_id2cumu_size has an extra fake chromosome (0) compared to chr_id2size
		data_ls = get_chr_id2cumu_size(self.chr_id2size)
		self.chr_id2cumu_size, self.chr_gap, self.chr_id_ls = data_ls[:3]
		#2008-02-04 update the info related to chromosome , position in toolbar
		if self.debug:
			print 'self.chr_id2size', self.chr_id2size
			for chr_id in self.chr_id2size:
				print type(chr_id)
			print 'self.chr_id2cumu_size', self.chr_id2cumu_size
			for chr_id in self.chr_id2cumu_size:
				print type(chr_id)
			print 'self.chr_id_ls', self.chr_id_ls
			for chr_id in self.chr_id_ls:
				print type(chr_id)
		self.toolbar.update_chr_info(self.chr_id2size, self.chr_id2cumu_size, self.chr_gap, self.chr_id_ls)
	
	def getXposition(self, chr, pos):
		chr = str(chr)
		if getattr(self, 'chr_id2cumu_size', None) is None:
			self.load_data(self.mysql_curs, self.postgres_curs)
		this_chr_starting_pos_on_plot = self.chr_id2cumu_size[chr]-self.chr_id2size[chr]-self.chr_gap
		x = this_chr_starting_pos_on_plot + pos
		return x
	
	def plot(self, ax, canvas, genome_wide_result, draw_line_as_point=True):
		"""
		2008-05-28
			input is genome_wide_result
			chr_id2cumu_size, chr_id2size, chr_gap hidden from arguments
		2008-02-04
			chromosome is converted to str type
		2008-02-01
			draw the p-value, snp position, and chromosome boundary
		"""
		sys.stderr.write("Plotting %s ..."%genome_wide_result.name)
		#ax.clear()
		genome_wide_result_id = id(genome_wide_result)
		x_ls = []
		y_ls = []
		for data_obj in genome_wide_result.data_obj_ls:
			y_pos = genome_wide_result.base_value - genome_wide_result.min_value + data_obj.value 
			x_pos = self.getXposition(data_obj.chromosome, data_obj.position)
			if data_obj.stop_position is not None and draw_line_as_point==False:	#bigger than 100k, then a line
				x_stop_pos = self.getXposition(data_obj.chromosome, data_obj.stop_position)
				x_ls.append([(x_pos, y_pos), (x_stop_pos, y_pos)])
				#artist_obj = Line2D([x_pos, y_pos], [x_stop_pos, y_pos], picker=True)
			else:
				if draw_line_as_point and data_obj.stop_position is not None:
					x_stop_pos = self.getXposition(data_obj.chromosome, data_obj.stop_position)
					x_pos = (x_pos+x_stop_pos)/2.0
				x_ls.append(x_pos)
				y_ls.append(y_pos)
					#artist_obj = Circle((x_pos, y_pos), picker=True)
		
		if len(y_ls)>0:
			artist_obj = ax.scatter(x_ls, y_ls, s=10, faceted=False, picker=True)
		else:
			artist_obj = LineCollection(x_ls, picker=True)
			ax.add_artist(artist_obj)
		artist_obj_id = id(artist_obj)
		self.artist_obj_id2data_obj_key[artist_obj_id] = [genome_wide_result_id, None]
		
		y_base_value = genome_wide_result.base_value
		y_top_value = genome_wide_result.base_value + genome_wide_result.max_value - genome_wide_result.min_value
		if self.debug:
			print "y_base_value", y_base_value
			print 'y_top_value', y_top_value
		self.yticks.append(y_base_value)
		self.yticks.append(y_top_value)
		ax.set_yticks(self.yticks)
		
		self.yticklabels.append('%s %.2f'%(genome_wide_result.name, genome_wide_result.min_value))
		self.yticklabels.append('%s %.2f'%(genome_wide_result.name, genome_wide_result.max_value))
		ax.set_yticklabels(self.yticklabels)
		
		"""
		ax.add_artist(g_artist)
		artist_obj_id = id(g_artist)
				self.artist_obj_id2artist_gene_id_ls[artist_obj_id] = [g_artist, gene_id]
				self.gene_id2artist_object_id[gene_id] = artist_obj_id
				
				x_ls = []
		y_ls = []
		max_pvalue = 0
		for i in range(len(snp_pos_ls)):
			chr, pos = snp_pos_ls[i]
			chr = str(chr)
			this_chr_starting_pos_on_plot = chr_id2cumu_size[chr]-chr_id2size[chr]-chr_gap
			x = this_chr_starting_pos_on_plot + pos
			x_ls.append(x)
			pvalue = pvalue_ls[i]
			if pvalue > max_pvalue:
				max_pvalue = pvalue
			y_ls.append(pvalue)
		ax.plot(x_ls, y_ls, '.', picker=3)	#3 points tolerance
		"""
		#draw the chromosome boundary
		for chr_id, cumu_size in self.chr_id2cumu_size.iteritems():
			ax.vlines(cumu_size, y_base_value, y_top_value, color='k')
		canvas.draw()
		sys.stderr.write("Done.\n")
	
	def respond2GeneObjPicker(self, event, artist_obj_id, artist_obj_id2artist_gene_id_ls, gene_annotation,\
							ax=None, draw_gene_symbol_when_clicked=False, canvas_matplotlib=None):
		"""
		2008-12-17
			a common function specifying how to respond to picker events of 
				both gene models in axe_gene_model and vertical gene spans in ax
			called by on_canvas_pick()
		"""
		gene_id = artist_obj_id2artist_gene_id_ls[artist_obj_id][1]
		gene_model = gene_annotation.gene_id2model[gene_id]
		
		if len(gene_model.gene_commentaries)==0:
			gene_commentary = gene_model	#fake one here
		else:
			gene_commentary = gene_model.gene_commentaries[0]
		
		protein_label = getattr(gene_commentary, 'protein_label', None)
		if not protein_label:
			protein_label = getattr(gene_commentary, 'label', '')
		
		protein_comment = getattr(gene_commentary, 'protein_comment', None)
		if not protein_comment:
			protein_comment = getattr(gene_commentary, 'comment', '')
		
		if getattr(gene_commentary, 'protein_label', None) is not None:	#true gene_commentary is available
			type_of_gene = getattr(gene_model, 'type_of_gene', '')
		else:	#it doesn't have protein, get gene_commentary_type
			type_of_gene = getattr(gene_commentary, 'gene_commentary_type', '')
		
		print '%s (gene id=%s) type_of_gene: %s. chromosome: %s. start: %s. stop: %s. strand: %s.'%\
				(gene_model.gene_symbol, gene_id, type_of_gene, gene_model.chromosome, \
				gene_model.start, gene_model.stop, gene_model.strand)
		print '\t protein_label: %s.'%protein_label
		print '\t protein_comment: %s.'%protein_comment
		
		if draw_gene_symbol_when_clicked:
			if ax:
				ax.text(event.mouseevent.xdata, event.mouseevent.ydata, gene_model.gene_symbol, size=8)
			if canvas_matplotlib:
				canvas_matplotlib.draw()
	
	def on_canvas_pick(self, event):
		"""
		2008-11-12
			display more information (maf, genotype_var_perc, comment) of data_obj (SNP) if they exist
		2008-05-28
			pick from collection
		2008-01-31 copied from examples/pick_event_demo.py from matplotlib source code
		"""
		if self.debug:
			print dir(event)
			print event.artist
			print dir(event.artist)
			print type(event.artist)
		"""
		if isinstance(event.artist, Line2D):
			thisline = event.artist
			xdata = thisline.get_xdata()
			ydata = thisline.get_ydata()
			ind = event.ind
			if self.debug:
				print "indices:", ind
				print 'onpick1 line:', zip(numpy.take(xdata, ind), numpy.take(ydata, ind))
			for i in ind:
				print "snp chromosome: %s, position: %s, pvalue: %s"%(self.snp_pos_ls[i][0], self.snp_pos_ls[i][1], self.pvalue_ls[i])
		"""
		if isinstance(event.artist, ExonIntronCollection):	#ExonIntronCollection is also a kind of Collection
			artist_obj_id = id(event.artist)
			if artist_obj_id in self.artist_obj_id2artist_gene_id_ls:
				self.respond2GeneObjPicker(event, artist_obj_id, self.artist_obj_id2artist_gene_id_ls, self.gene_annotation,\
							ax=self.axe_gene_model, draw_gene_symbol_when_clicked=self.draw_gene_symbol_when_clicked, \
							canvas_matplotlib=self.canvas_matplotlib)
			else:
				sys.stderr.write("%s not in artist_obj_id2artist_gene_id_ls.\n"%(artist_obj_id))
		elif isinstance(event.artist, Collection) or isinstance(event.artist, LineCollection):	#
			artist_obj_id = id(event.artist)
			if artist_obj_id in self.artist_obj_id2data_obj_key:
				genome_wide_result_id, data_obj_id = self.artist_obj_id2data_obj_key[artist_obj_id]
				genome_wide_result = self.genome_wide_results.get_genome_wide_result_by_obj_id(genome_wide_result_id)
				for obj_index in event.ind:
					if isinstance(obj_index, tuple) or isinstance(obj_index, list):
						obj_index = obj_index[0]
					data_obj = genome_wide_result.get_data_obj_by_obj_index(obj_index)
					output_str = "genome result: %s, chromosome: %s, position: %s, "%(genome_wide_result.name, data_obj.chromosome, data_obj.position)
					if data_obj.stop_position is not None:
						output_str += "stop position: %s, "%(data_obj.stop_position)
					output_str += '\n'
					output_str += "\tscore: %s\n"%(data_obj.value)
					if data_obj.maf:
						output_str += "\tmaf: %s\n"%(data_obj.maf)
					if data_obj.genotype_var_perc:
						output_str += "\tgenotype_var_perc: %s\n"%(data_obj.genotype_var_perc)
					if data_obj.comment:
						output_str += "\tcomment: %s\n"%(data_obj.comment)
					print output_str
			else:
				sys.stderr.write("%s not in artist_obj_id2data_obj_key.\n"%(artist_obj_id))
		elif isinstance(event.artist, Polygon):	#ExonIntronCollection is also a kind of Collection
			artist_obj_id = id(event.artist)
			if artist_obj_id in self.artist_obj_id2artist_gene_id_ls:
				self.respond2GeneObjPicker(event, artist_obj_id, self.artist_obj_id2artist_gene_id_ls, self.gene_annotation,\
										ax=self.ax, draw_gene_symbol_when_clicked=self.draw_gene_symbol_when_clicked, \
										canvas_matplotlib=self.canvas_matplotlib)
			else:
				sys.stderr.write("%s not in artist_obj_id2artist_gene_id_ls.\n"%(artist_obj_id))
		elif isinstance(event.artist, Rectangle):
			patch = event.artist
			print 'onpick1 patch:', patch.get_verts()
		elif isinstance(event.artist, Text):
			text = event.artist
			print 'onpick1 text:', text.get_text()
		
	
	def on_imagemenuitem_quit_activate(self, data=None):
		"""
		2008-02-01
			program quits
		"""
		gtk.main_quit()
	
	def on_imagemenuitem_open_activate(self, event, data=None):
		self.filechooserdialog1.show_all()
	
	def on_imagemenuitem_db_connect_activate(self, event, data=None):
		self.dialog_db_connect.show_all()
	
	def on_button_filechooser_ok_clicked(self, widget, data=None):
		"""
		2008-12-16
			allow gwr name to be specified
			add function to get gwr from db based on call_method_id, analysis_method_id, phenotype_method_id
		2008-10-12
			add checkbutton_draw_line_as_point
			add checkbutton_4th_col_stop_pos
		2008-08-03
			restrict the data by (chromosome, start, stop)
		2008-05-31
			add check button to handle log10 transformation
		2008-05-28
			use GenomeWideResult and etc
		2008-02-14
			set the window title by the input filename
		"""
		input_fname = self.filechooserdialog1.get_filename()
		self.filechooserdialog1.hide()
		if not self.mysql_conn or not self.mysql_curs:
			self.db_connect()
		self.app1.set_title("Genome Browser: %s"%input_fname)
		
		checkbutton_log10_transformation = self.xml.get_widget("checkbutton_log10_transformation")
		if checkbutton_log10_transformation.get_active():
			do_log10_transformation = True
		else:
			do_log10_transformation = False
		
		if self.entry_min_value_cutoff.get_text():
			min_value_cutoff = float(self.entry_min_value_cutoff.get_text())
		else:
			min_value_cutoff = None
		#2008-08-03
		pdata = PassingData()
		entry_chromosome = self.xml.get_widget("entry_chromosome")
		if entry_chromosome.get_text():
			pdata.chromosome = int(entry_chromosome.get_text())
		entry_start = self.xml.get_widget("entry_start")
		if entry_start.get_text():
			pdata.start = int(entry_start.get_text())
		entry_stop = self.xml.get_widget("entry_stop")
		if entry_stop.get_text():
			pdata.stop = int(entry_stop.get_text())
		
		checkbutton_4th_col_stop_pos = self.xml.get_widget("checkbutton_4th_col_stop_pos")
		if checkbutton_4th_col_stop_pos.get_active():
			pdata.is_4th_col_stop_pos = True
		else:
			pdata.is_4th_col_stop_pos = False
		
		checkbutton_draw_line_as_point = self.xml.get_widget("checkbutton_draw_line_as_point")
		if checkbutton_draw_line_as_point.get_active():
			draw_line_as_point= True
		else:
			draw_line_as_point = False
		
		entry_gwr_name = self.xml.get_widget("entry_gwr_name")
		if entry_gwr_name.get_text():
			pdata.gwr_name = entry_gwr_name.get_text()
		else:
			pdata.gwr_name = None
		
		entry_call_method_id = self.xml.get_widget("entry_call_method_id")
		call_method_id = entry_call_method_id.get_text()
		entry_analysis_method_id = self.xml.get_widget("entry_analysis_method_id")
		analysis_method_id = entry_analysis_method_id.get_text()
		entry_phenotype_method_id = self.xml.get_widget("entry_phenotype_method_id")
		phenotype_method_id = entry_phenotype_method_id.get_text()
		
		if call_method_id and analysis_method_id and phenotype_method_id:
			call_method_id = int(call_method_id)
			analysis_method_id = int(analysis_method_id)
			phenotype_method_id = int(phenotype_method_id)
			rows = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).filter_by(analysis_method_id=analysis_method_id).\
					filter_by(phenotype_method_id=phenotype_method_id).filter_by(results_method_type_id=1)
			if rows.count()==1:
				rm = rows.first()
			elif rows.count()==0:
				sys.stderr.write("No result fetched from db based on call_method_id=%s, analysis_method_id=%s, phenotype_method_id=%s.\n"%\
								(call_method_id, analysis_method_id, phenotype_method_id))
				rm = None
			else:
				sys.stderr.write("First result out of %s results fetched from db based on call_method_id=%s, analysis_method_id=%s, phenotype_method_id=%s.\n"%\
								(rows.count(), call_method_id, analysis_method_id, phenotype_method_id))
				rm = rows.first()
			if rm:
				input_fname = rm.filename
				pdata.gwr_name = '%s_%s_%s'%(rm.analysis_method.short_name, rm.phenotype_method_id, rm.phenotype_method.short_name)
			
		
		genome_wide_result = getGenomeWideResultFromFile(input_fname, min_value_cutoff, do_log10_transformation, pdata)
		if len(genome_wide_result.data_obj_ls)>0:
			self.genome_wide_results.add_genome_wide_result(genome_wide_result)
			#self.load_data(input_fname, self.mysql_curs, self.postgres_curs)
			self.plot(self.ax, self.canvas_matplotlib, self.genome_wide_results.genome_wide_result_ls[-1], draw_line_as_point=draw_line_as_point)
		else:
			sys.stderr.write("No data in %s under min_value_cutoff=%s. Maybe min_value_cutoff is too high.\n"%(input_fname, min_value_cutoff))
	
	def on_button_filechooser_cancel_clicked(self, widget, data=None):
		self.filechooserdialog1.hide()
	
	def on_button_dialog_db_connect_cancel_clicked(self, widget, data=None):
		self.dialog_db_connect.hide()
	
	def on_button_dialog_db_connect_clicked(self, widget, data=None):
		self.dialog_db_connect.hide()
		self.db_connect()
	
	def db_connect(self):
		"""
		2008-12-16
			add gene_annotation_picklef
		2008-02-01
			read the data in dialog_db_connect and establish the connections to two databases
		"""
		sys.stderr.write("Database Connecting ...")
		hostname = self.entry_mysql_hostname.get_text()
		dbname = self.entry_mysql_dbname.get_text()
		db_user = self.xml.get_widget("entry_db_user").get_text()
		db_passwd = self.xml.get_widget("entry_db_passwd").get_text()
		
		import MySQLdb
		try:
			self.mysql_conn = MySQLdb.connect(db=dbname,host=hostname)
			self.mysql_curs = self.mysql_conn.cursor()
			self.db = Stock_250kDB.Stock_250kDB(drivername='mysql', username=db_user,
					   password=db_passwd, hostname=hostname, database=dbname)
			self.db.setup(create_tables=False)
			self.session = self.db.session
		except:
			sys.stderr.write('DB connection error: %s\n'%repr(sys.exc_info()))
			traceback.print_exc()
		
		hostname = self.entry_postgres_hostname.get_text()
		dbname = self.entry_postgres_dbname.get_text()
		schema = self.entry_postgres_schema.get_text()
		
		if not self.gene_annotation:
			gene_annotation_picklef = self.entry_gene_annotation_picklef.get_text()
			self.gene_annotation = DrawSNPRegion.dealWithGeneAnnotation(gene_annotation_picklef)
		
		#from annot.bin.codense.common import db_connect			#2008-12-16 don't need postgres conn anymore
		#self.postgres_conn, self.postgres_curs = db_connect(hostname, dbname, schema)
		
		sys.stderr.write("Done.\n")
	
	def get_gene_id2model(cls, curs, entrezgene_mapping_table='genome.entrezgene_mapping', \
						annot_assembly_table = 'genome.annot_assembly', gene_table='genome.gene', \
						gene2go_table='genome.gene2go', tax_id=3702):
		"""
		2008-09-24
			turn gene_id into integer
		2008-08-03
			schema where tables about genes are from is renamed from 'sequence' to 'genome'
		2008-02-02
			get all the necessary info for genes.
			watch, chromosome here is varchar type (because of chromosome X, Y etc)
		"""
		sys.stderr.write("Getting gene_id2model and chr_id2gene_id_ls...")
		from annot.bin.codense.common import pg_1d_array2python_ls
		gene_id2model = {}
		chr_id2gene_id_ls = {}
		curs.execute("DECLARE gene_crs CURSOR FOR select e.gene_id, a.chromosome, e.start, e.stop, e.mrna_start, e.mrna_stop, e.cds_start, e.cds_stop, e.strand, g.gene_symbol, g.description, g.type_of_gene \
					from %s e, %s a, %s g where e.gene_id=g.gene_id and e.genomic_gi=a.gi and e.tax_id=%s order by chromosome, start, stop"%\
					(entrezgene_mapping_table, annot_assembly_table, gene_table, tax_id))
		curs.execute("fetch 5000 from gene_crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				#gene_id is integer. chromosome is varchar.
				gene_id, chromosome, start, stop, mrna_start, mrna_stop, cds_start, cds_stop, strand, symbol, description, type_of_gene = row
				gene_id = int(gene_id)	#2008-09-24
				if cds_start and cds_stop:
					if type(cds_start)!=list:
						cds_start = pg_1d_array2python_ls(cds_start, int)
						cds_stop = pg_1d_array2python_ls(cds_stop, int)
				else:
					cds_start = cds_stop = None
				
				if mrna_start and mrna_stop:
					if type(mrna_stop)!=list:
						mrna_start = pg_1d_array2python_ls(mrna_start, int)
						mrna_stop = pg_1d_array2python_ls(mrna_stop, int)
				else:
					mrna_start = mrna_stop = None
				
				if chromosome not in chr_id2gene_id_ls:
					chr_id2gene_id_ls[chromosome] = []
				chr_id2gene_id_ls[chromosome].append(gene_id)
				if gene_id not in gene_id2model:
					gene_id2model[gene_id] = GeneModel(gene_id, chromosome, symbol, description, type_of_gene, \
														start, stop, mrna_start, mrna_stop, cds_start, cds_stop, strand)
			curs.execute("fetch 5000 from gene_crs")
			rows = curs.fetchall()
		curs.execute("close gene_crs")
		sys.stderr.write("Done.\n")
		return gene_id2model, chr_id2gene_id_ls
	
	get_gene_id2model = classmethod(get_gene_id2model)
	
	def plot_one_gene(self, ax, gene_id, gene_id2model, chr_id2cumu_size, chr_id2size, chr_gap, y_value=1, gene_width=1.0):
		"""
		2008-12-16
			defunct. DrawSNPRegion.drawGeneModel() is used in on_button_draw_annotation_clicked()
		2008-02-02
			draw a single gene on the canvas, 
		"""
		gene_model = gene_id2model.get(gene_id)
		if gene_model:
			c_start_ls = None
			c_end_ls = None
			if gene_model.cds_start!=None and gene_model.cds_stop!=None:
				c_start_ls = gene_model.cds_start
				c_end_ls = gene_model.cds_stop
			elif gene_model.mrna_start!=None and gene_model.mrna_stop!=None:
				c_start_ls = gene_model.mrna_start
				c_end_ls = gene_model.mrna_stop
			elif gene_model.start!=None and gene_model.stop!=None:
				c_start_ls = [gene_model.start]
				c_end_ls = [gene_model.stop]
			if c_start_ls and c_end_ls:
				chromosome = gene_model.chromosome
				this_chr_starting_pos_on_plot = chr_id2cumu_size[chromosome]-chr_id2size[chromosome]-chr_gap
				if gene_model.strand=="1":
					g_artist = Gene(c_start_ls, c_end_ls, y=y_value, x_offset=this_chr_starting_pos_on_plot, width=gene_width, alpha=0.3, facecolor='r', picker=True)
				elif gene_model.strand=="-1":	#to draw opposite strand, 1st is to order c_start_ls and c_end_ls in descending order. 2nd is to swap c_start_ls and c_end_ls.
					#c_start_ls.reverse()	#2008-02-04 it's already in descending order in db.
					#c_end_ls.reverse()	#2008-02-04 it's already in descending order in db.
					g_artist = Gene(c_end_ls, c_start_ls, y=y_value, x_offset=this_chr_starting_pos_on_plot, width=gene_width, alpha=0.3, facecolor='r', picker=True)
				else:	#no arrow
					g_artist = Gene(c_start_ls, c_end_ls, y=y_value, is_arrow=False, x_offset=this_chr_starting_pos_on_plot, width=gene_width, alpha=0.3, facecolor='r', picker=True)
				ax.add_artist(g_artist)
				artist_obj_id = id(g_artist)
				self.artist_obj_id2artist_gene_id_ls[artist_obj_id] = [g_artist, gene_id]
				self.gene_id2artist_object_id[gene_id] = artist_obj_id
	
	def on_button_draw_annotation_clicked(self, widget, data=None):
		"""
		2008-12-16
			use DrawSNPRegion.drawGeneModel() to draw gene models
		2008-02-02
		"""
		if not self.chr_id2size:
			sys.stderr.write("No genome-wide pvalue plot has been drawn yet. Do it first!\n")
			return
		#if not self.gene_id2model:
		#	self.gene_id2model, self.chr_id2gene_id_ls = self.get_gene_id2model(self.postgres_curs, tax_id=3702)
		if not self.gene_annotation:
			self.gene_annotation = DrawSNPRegion.dealWithGeneAnnotation(self.entry_gene_annotation_picklef.get_text())
		
		xlim = self.axe_gene_model.get_xlim()
		left_chr, left_pos = get_chr_pos_from_x_axis_pos(xlim[0], self.chr_gap, self.chr_id2cumu_size, self.chr_id_ls)
		right_chr, right_pos = get_chr_pos_from_x_axis_pos(xlim[1], self.chr_gap, self.chr_id2cumu_size, self.chr_id_ls)
		
		#fake a snps_within_this_region for drawGeneModel()
		snps_within_this_region = PassingData(chr_pos_ls=[[left_chr, left_pos],[right_chr, right_pos]])
		base_y_value = 1
		gene_width = 0.8
		gene_position_cycle = 5
		
		return_data = DrawSNPRegion.drawGeneModel(self.axe_gene_model, snps_within_this_region, self.gene_annotation, candidate_gene_set=None,\
								gene_width=gene_width, gene_position_cycle=gene_position_cycle, base_y_value=base_y_value, \
								gene_box_text_gap=20, label_gene=0, rotate_xy=False,\
								chr_id2cumu_size=self.chr_id2cumu_size, chr_id2size=self.chr_id2size, chr_gap=self.chr_gap,\
								artist_obj_id2artist_gene_id_ls=self.artist_obj_id2artist_gene_id_ls, \
								gene_id2artist_object_id=self.gene_id2artist_object_id, drawGeneOnTheBoundary=False)
					#set drawGeneOnTheBoundary to False because later adding text to these genes would corrupt the running program.
		self.axe_gene_model.set_ylim([base_y_value-gene_width, gene_position_cycle+gene_width*2])
		
		"""
		for gene_id in self.chr_id2gene_id_ls[left_chr]:
			gene_model = self.gene_id2model[gene_id]
			if gene_model.start!=None and gene_model.stop!=None and gene_model.stop>left_pos and gene_id not in self.gene_id2artist_object_id:
				if left_chr==right_chr:	#same chromosome
					if gene_model.start>right_pos:	#totally out of range, skip it
						continue
				y_value = len(self.gene_id2artist_object_id)%4	#cycling through the y position to avoid clogging
				self.plot_one_gene(self.ax, gene_id, self.gene_id2model, self.chr_id2cumu_size, self.chr_id2size, self.chr_gap, y_value=-1-y_value, gene_width=self.gene_width)
		if left_chr!=right_chr:
			for gene_id in self.chr_id2gene_id_ls[right_chr]:
				gene_model = self.gene_id2model[gene_id]
				if gene_model.start!=None and gene_model.stop!=None and gene_model.start<right_pos and gene_id not in self.gene_id2artist_object_id:
					y_value = len(self.gene_id2artist_object_id)%4	#cycling through the y position to avoid clogging
					self.plot_one_gene(self.ax, gene_id, self.gene_id2model, self.chr_id2cumu_size, self.chr_id2size, self.chr_gap, y_value=-1-y_value, gene_width=self.gene_width)
		"""
		self.canvas_matplotlib.draw()
	
	def on_imagemenuitem_preferences_activate(self, event, data=None):
		"""
		2008-02-04
		"""
		self.dialog_preferences.show_all()
	
	def on_button_dialog_preferences_ok_clicked(self, widget, data=None):
		"""
		2008-02-04
			change some preferences
		"""
		self.dialog_preferences.hide()
		if self.checkbutton_debug.get_active():
			self.debug = 1
		else:
			self.debug = 0
		if self.checkbutton_stderr.get_active():
			sys.stderr = self.dummy_err
		else:
			sys.stderr = sys.__stderr__
		if self.checkbutton_stdout.get_active():
			sys.stdout = self.dummy_out
		else:
			sys.stdout = sys.__stdout__
		if self.checkbutton_draw_gene_symbol.get_active():
			self.draw_gene_symbol_when_clicked = 1
		else:
			self.draw_gene_symbol_when_clicked = 0
		self.gene_width = float(self.entry_gene_width.get_text())
	
	def on_button_dialog_preferences_cancel_clicked(self, widget, data=None):
		"""
		2008-02-04
			don't change any preferences
		"""
		self.dialog_preferences.hide()
	
	def on_imagemenuitem_about_activate(self, widget):
		"""
		2008-02-04
		"""
		self.aboutdialog1.show_all()
	
	def on_imagemenuitem_cleanup_output_activate(self, widget):
		"""
		2008-02-04
			clean up output buffer
		"""
		self.textbuffer_output.set_text('')
	
	def on_checkbutton_debug_toggled(self, widget):
		"""
		2008-05-28
		"""
		if self.checkbutton_debug.get_active():
			self.debug = 1
		else:
			self.debug = 0
	
	def on_checkbutton_stdout_toggled(self, widget):
		"""
		2008-05-28
		"""
		if self.checkbutton_stdout.get_active():
			sys.stdout = self.dummy_out
		else:
			sys.stdout = sys.__stdout__
	
	def on_checkbutton_stderr_toggled(self, widget):
		"""
		2008-05-28
		"""
		if self.checkbutton_stderr.get_active():
			sys.stderr = self.dummy_err
		else:
			sys.stderr = sys.__stderr__
	
	def on_checkbutton_draw_gene_symbol_toggled(self, widget):
		"""
		2008-05-28
		"""
		if self.checkbutton_draw_gene_symbol.get_active():
			self.draw_gene_symbol_when_clicked = 1
		else:
			self.draw_gene_symbol_when_clicked = 0
	
	def on_entry_gene_width_changed(self, widget):
		"""
		2008-05-28
		"""
		self.gene_width = float(self.entry_gene_width.get_text())
	
	def on_filechooserbutton_gene_annot_file_set(self, widget):
		"""
		2008-12-16
		"""
		self.entry_gene_annotation_picklef.set_text(self.filechooserbutton_gene_annot.get_filename())
	
	def on_button_draw_gene_list_bars_clicked(self, widget):
		"""
		2008-12-16
			draw vertical spans to denote the locations of genes from a candidate list
		"""
		if self.db is None:
			self.db_connect()
		if not self.chr_id2size:
			sys.stderr.write("No genome-wide pvalue plot has been drawn yet. Do it first!\n")
			return
		entry_gene_list_id = self.xml.get_widget("entry_gene_list_id")
		list_type_id = entry_gene_list_id.get_text()
		comboboxentry_bar_color = self.xml.get_widget("comboboxentry_bar_color")
		bar_color = comboboxentry_bar_color.get_active_text()
		if not bar_color:	#default is black
			bar_color = 'k'
		if list_type_id:
			list_type_id = int(list_type_id)
			self.candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
			for gene_id in self.candidate_gene_set:
				gene_model = self.gene_annotation.gene_id2model[gene_id]
				if gene_id in self.gene_id2vspan_obj_id:
					artist_obj_id = self.gene_id2vspan_obj_id[gene_id]
					artist = self.artist_obj_id2artist_gene_id_ls[artist_obj_id][0]
					if artist.get_edgecolor()!=bar_color:
						artist.set_edgecolor(bar_color)
					if artist.get_facecolor()!=bar_color:
						artist.set_facecolor(bar_color)
					#artist.remove()
				else:
					this_chr_starting_pos_on_plot = self.chr_id2cumu_size[gene_model.chromosome]-\
							self.chr_id2size[gene_model.chromosome]-self.chr_gap
					xmin = this_chr_starting_pos_on_plot + gene_model.start
					xmax = this_chr_starting_pos_on_plot + gene_model.stop
					artist = self.ax.axvspan(xmin, xmax, edgecolor=bar_color, facecolor=bar_color, alpha=0.3, picker=6)
					artist_obj_id = id(artist)
					self.artist_obj_id2artist_gene_id_ls[artist_obj_id] = [artist, gene_id]
					self.gene_id2vspan_obj_id[gene_id] = artist_obj_id
			self.canvas_matplotlib.draw()
	
	def on_button_adjust_gene_axis_clicked(self, widget):
		"""
		2008-12-19
			sometimes after zoom-in/out, axe_gene_model loses track of its y-range and the gene models in it float into ax.
			this function would bring the y-range of axe_gene_model into normal range.
		"""
		base_y_value = 1
		gene_width = 0.8
		gene_position_cycle = 5
		self.axe_gene_model.set_ylim([base_y_value-gene_width, gene_position_cycle+gene_width*2])
		self.canvas_matplotlib.draw()
	
if __name__ == '__main__':
	prog = gnome.program_init('GenomeBrowser', '0.1')
	instance = GenomeBrowser()
	gtk.main()
