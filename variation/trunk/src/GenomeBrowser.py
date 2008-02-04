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
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import pymodule.gnome as yh_gnome

import numpy
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.text import Text

from pymodule.yh_matplotlib_artists import Gene

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

class GenomeBrowser:
	def __init__(self):
		"""
		2008-01-30
		"""
		program_path = os.path.dirname(sys.argv[0])
		xml = gtk.glade.XML(os.path.join(program_path, 'GenomeBrowser.glade'))
		xml.signal_autoconnect(self)
		self.app1 = xml.get_widget("app1")
		self.app1.connect("delete_event", gtk.main_quit)
		self.app1.set_default_size(1200, 800)
		
		self.vbox_matplotlib = xml.get_widget('vbox_matplotlib')
		
		# matplotlib canvas
		fig = Figure(figsize=(8,8))
		self.canvas_matplotlib = FigureCanvas(fig)  # a gtk.DrawingArea
		self.canvas_matplotlib.set_size_request(800,600)
		self.canvas_matplotlib.mpl_connect('pick_event', self.on_canvas_pick)
		self.vbox_matplotlib.pack_start(self.canvas_matplotlib)
		
		#matplotlib axes
		self.ax = fig.add_subplot(111)
		
		# matplotlib toolbar
		toolbar = NavigationToolbar(self.canvas_matplotlib, self.app1)
		self.vbox_matplotlib.pack_start(toolbar, False, False)
		
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
		
		self.aboutdialog1 = xml.get_widget("aboutdialog1")
		self.aboutdialog1.connect("delete_event", yh_gnome.subwindow_hide)
		
		self.mysql_conn = self.mysql_curs = self.postgres_conn = self.postgres_curs = None
		
		self.chr_id2size = None
		self.chr_id2cumu_size = None
		
		self.gene_id2artist_object_id = {}
		self.chr_id2gene_id_ls = {}	#chr_id here is str type (db is varchar type)
		self.gene_id2model = {}
		self.artist_object_id2artist_gene_id_ls = {}
		
		self.gene_width = 1.0
		
		self.debug = 1
		
	def load_data(self, input_fname, mysql_curs, postgres_curs):
		"""
		2008-02-01
			read the input data
			
			chr_id2cumu_size has an extra fake chromosome (0) compared to chr_id2size
			
			chr_id is all changed into str type
		"""
		sys.stderr.write("Read in data from %s ... "%input_fname)
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		self.snp_pos_ls = []
		self.pvalue_ls = []
		for row in reader:
			chr, pos, pvalue = row
			self.snp_pos_ls.append([int(chr), int(pos)])
			self.pvalue_ls.append(float(pvalue))
		del reader
		sys.stderr.write("Done.\n")
		
		from variation.src.common import get_chr_id2size, get_chr_id2cumu_size
		chr_id_int2size = get_chr_id2size(mysql_curs)
		self.chr_id2size = {}	#change the type of chr_id into string type
		for chr_id_int, size in chr_id_int2size.iteritems():
			self.chr_id2size[str(chr_id_int)] = size
		#chr_id2cumu_size has an extra fake chromosome (0) compared to chr_id2size
		data_ls = get_chr_id2cumu_size(self.chr_id2size)
		self.chr_id2cumu_size, self.chr_gap, self.chr_id_ls = data_ls[:3]
		
	def plot(self, ax, canvas, snp_pos_ls, pvalue_ls, chr_id2cumu_size, chr_id2size, chr_gap):
		"""
		2008-02-04
			chromosome is converted to str type
		2008-02-01
			draw the p-value, snp position, and chromosome boundary
		"""
		sys.stderr.write("Plotting ...")
		ax.clear()
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
		
		#draw the chromosome boundary
		for chr_id, cumu_size in chr_id2cumu_size.iteritems():
			ax.plot([cumu_size, cumu_size], [-2, max_pvalue], c='k')
		canvas.draw()
		sys.stderr.write("Done.\n")
	
	def on_canvas_pick(self, event):
		"""
		2008-01-31 copied from examples/pick_event_demo.py from matplotlib source code
		"""
		if self.debug:
			print dir(event)
			print event.artist
			print dir(event.artist)
			print type(event.artist)
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
			print 'gene symbol: %s. description: %s. type_of_gene: %s. chromosome: %s. start: %s. stop: %s. strand: %s.'%\
				(gene_model.symbol, gene_model.description, gene_model.type_of_gene, gene_model.chromosome, gene_model.start, gene_model.stop, gene_model.strand)
	
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
		input_fname = self.filechooserdialog1.get_filename()
		self.filechooserdialog1.hide()
		if not self.mysql_conn or not self.mysql_curs or not self.postgres_conn or not self.postgres_curs:
			self.db_connect()
		self.load_data(input_fname, self.mysql_curs, self.postgres_curs)
		self.plot(self.ax, self.canvas_matplotlib, self.snp_pos_ls, self.pvalue_ls, self.chr_id2cumu_size, self.chr_id2size, self.chr_gap)
	
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
	
	def get_gene_id2model(self, curs, entrezgene_mapping_table='sequence.entrezgene_mapping', \
						annot_assembly_table = 'annot_assembly', gene_table='gene.gene', \
						gene2go_table='gene.gene2go', tax_id=3702):
		"""
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
	
	def get_chr_pos_from_x_axis_pos(self, x_axis_pos, chr_gap, chr_id2cumu_size, chr_id_ls):
		"""
		2008-02-03
			get chromosome, position from the x axis position
			chr_id_ls is the sorted version of the keys of chr_id2cumu_size
		"""
		chr_id_chosen = 0
		position = -1
		for i in range(1, len(chr_id_ls)):	#the 1st in chr_id_ls is fake chromosome 0
			prev_chr_id = chr_id_ls[i-1]
			chr_id = chr_id_ls[i]
			if chr_id2cumu_size[chr_id]>=x_axis_pos:
				chr_id_chosen = chr_id
				position = x_axis_pos - chr_id2cumu_size[prev_chr_id]
				break
		return chr_id_chosen, position
	
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
				elif gene_model.strand=="-1":
					c_start_ls.reverse()
					c_end_ls.reverse()
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
		left_chr, left_pos = self.get_chr_pos_from_x_axis_pos(xlim[0], self.chr_gap, self.chr_id2cumu_size, self.chr_id_ls)
		right_chr, right_pos = self.get_chr_pos_from_x_axis_pos(xlim[1], self.chr_gap, self.chr_id2cumu_size, self.chr_id_ls)
		
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
		self.gene_width = float(self.entry_gene_width.get_text())
	
	def on_button_dialog_preferences_cancel_clicked(self, widget, data=None):
		"""
		2008-02-04
			don't change any preferences
		"""
		self.dialog_preferences.hide()
	
	def on_imagemenuitem10_activate(self, widget):
		"""
		2008-02-04
		"""
		self.aboutdialog1.show_all()
	
prog = gnome.program_init('GenomeBrowser', '0.1')
instance = GenomeBrowser()
gtk.main()
