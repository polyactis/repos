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
from pymodule.gnome import NavigationToolbar2GTKAgg_chromosome as NavigationToolbar

import pymodule.gnome as yh_gnome

import numpy
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.text import Text
from matplotlib.collections import LineCollection, Collection

from pymodule.yh_matplotlib_artists import Gene
from variation.src.common import get_chr_pos_from_x_axis_pos
from pymodule.db import TableClass
from Results2DB_250k import Results2DB_250k
from pymodule import GenomeWideResults, GenomeWideResult, DataObject, getGenomeWideResultFromFile, PassingData

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
		self.ax = fig.add_subplot(111)
		
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
		
		self.dialog_preferences = xml.get_widget("dialog_preferences")
		self.dialog_preferences.connect("delete_event", yh_gnome.subwindow_hide)
		self.checkbutton_debug = xml.get_widget("checkbutton_debug")
		self.checkbutton_stdout = xml.get_widget("checkbutton_stdout")
		self.checkbutton_stderr = xml.get_widget("checkbutton_stderr")
		self.entry_gene_width = xml.get_widget("entry_gene_width")
		self.checkbutton_draw_gene_symbol = xml.get_widget("checkbutton_draw_gene_symbol")
		
		self.aboutdialog1 = xml.get_widget("aboutdialog1")
		self.aboutdialog1.connect("delete_event", yh_gnome.subwindow_hide)
		
		self.mysql_conn = self.mysql_curs = self.postgres_conn = self.postgres_curs = None
		
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
		self.artist_object_id2artist_gene_id_ls = {}
		
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
		artist_object_id = id(g_artist)
				self.artist_object_id2artist_gene_id_ls[artist_object_id] = [g_artist, gene_id]
				self.gene_id2artist_object_id[gene_id] = artist_object_id
				
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
	
	def on_canvas_pick(self, event):
		"""
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
		if isinstance(event.artist, Collection) or isinstance(event.artist, LineCollection):
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
					output_str += "value: %s"%(data_obj.value)
					print output_str
			else:
				sys.stderr.write("%s not in artist_obj_id2data_obj_key.\n"%(artist_obj_id))
		elif isinstance(event.artist, Rectangle):
			patch = event.artist
			print 'onpick1 patch:', patch.get_verts()
		elif isinstance(event.artist, Text):
			text = event.artist
			print 'onpick1 text:', text.get_text()
		elif isinstance(event.artist, Gene):
			artist_object_id = id(event.artist)
			gene_id = self.artist_object_id2artist_gene_id_ls[artist_object_id][1]
			gene_model = self.gene_id2model[gene_id]
			print 'gene id: %s. symbol: %s. description: %s. type_of_gene: %s. chromosome: %s. start: %s. stop: %s. strand: %s.'%\
				(gene_id, gene_model.symbol, gene_model.description, gene_model.type_of_gene, gene_model.chromosome, gene_model.start, gene_model.stop, gene_model.strand)
			if self.draw_gene_symbol_when_clicked:
				self.ax.text(event.mouseevent.xdata, event.mouseevent.ydata, gene_model.symbol, size=8)
				self.canvas_matplotlib.draw()
	
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
		if not self.mysql_conn or not self.mysql_curs or not self.postgres_conn or not self.postgres_curs:
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
		
		genome_wide_result = getGenomeWideResultFromFile(input_fname, min_value_cutoff, do_log10_transformation, pdata)
		if len(genome_wide_result.data_obj_ls)>0:
			self.genome_wide_results.add_genome_wide_result(genome_wide_result)
			#self.load_data(input_fname, self.mysql_curs, self.postgres_curs)
			self.plot(self.ax, self.canvas_matplotlib, self.genome_wide_results.genome_wide_result_ls[-1])
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
		2008-02-01
			read the data in dialog_db_connect and establish the connections to two databases
		"""
		sys.stderr.write("Database Connecting ...")
		hostname = self.entry_mysql_hostname.get_text()
		dbname = self.entry_mysql_dbname.get_text()
		import MySQLdb
		self.mysql_conn = MySQLdb.connect(db=dbname,host=hostname)
		self.mysql_curs = self.mysql_conn.cursor()
		
		hostname = self.entry_postgres_hostname.get_text()
		dbname = self.entry_postgres_dbname.get_text()
		schema = self.entry_postgres_schema.get_text()
		from annot.bin.codense.common import db_connect
		self.postgres_conn, self.postgres_curs = db_connect(hostname, dbname, schema)
		sys.stderr.write("Done.\n")
	
	def get_gene_id2model(self, curs, entrezgene_mapping_table='genome.entrezgene_mapping', \
						annot_assembly_table = 'genome.annot_assembly', gene_table='genome.gene', \
						gene2go_table='genome.gene2go', tax_id=3702):
		"""
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
				if cds_start and cds_stop:
					cds_start = pg_1d_array2python_ls(cds_start, int)
					cds_stop = pg_1d_array2python_ls(cds_stop, int)
				else:
					cds_start = cds_stop = None
				
				if mrna_start and mrna_stop:
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
	
	def plot_one_gene(self, ax, gene_id, gene_id2model, chr_id2cumu_size, chr_id2size, chr_gap, y_value=1, gene_width=1.0):
		"""
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
				artist_object_id = id(g_artist)
				self.artist_object_id2artist_gene_id_ls[artist_object_id] = [g_artist, gene_id]
				self.gene_id2artist_object_id[gene_id] = artist_object_id
	
	def on_button_draw_annotation_clicked(self, widget, data=None):
		"""
		2008-02-02
		"""
		if not self.chr_id2size:
			sys.stderr.write("No plot has been drawn yet. Open a file first!\n")
			return
		if not self.gene_id2model:
			self.gene_id2model, self.chr_id2gene_id_ls = self.get_gene_id2model(self.postgres_curs, tax_id=3702)
		from pymodule.yh_matplotlib_artists import Gene
		xlim = self.ax.get_xlim()
		left_chr, left_pos = get_chr_pos_from_x_axis_pos(xlim[0], self.chr_gap, self.chr_id2cumu_size, self.chr_id_ls)
		right_chr, right_pos = get_chr_pos_from_x_axis_pos(xlim[1], self.chr_gap, self.chr_id2cumu_size, self.chr_id_ls)
		
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

if __name__ == '__main__':
	prog = gnome.program_init('GenomeBrowser', '0.1')
	instance = GenomeBrowser()
	gtk.main()
