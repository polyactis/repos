#!/usr/bin/env python
"""
2008-01-17
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
import gtk, gtk.glade
from gtk import gdk
import gnome
import gnome.ui

import matplotlib
matplotlib.use('GTKAgg')  # or 'GTK'
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

from variation.src.QualityControl import QualityControl
from variation.src import Cmp250kVs2010, Cmp250kVs149SNP, Cmp149SNPVs2010, CmpTina2010VsMy2010In250kSNPs

class QCVisualizeII:
	def __init__(self):
		"""
		2008-01-17
		"""
		program_path = os.path.dirname(sys.argv[0])
		xml = gtk.glade.XML(os.path.join(program_path, 'QCVisualizeII.glade'))
		xml.signal_autoconnect(self)
		self.app1 = xml.get_widget("app1")
		self.app1.connect("delete_event", self.subwindow_hide)
		self.app1.set_default_size(800, 800)
		
		self.treeview_strains1 = xml.get_widget('treeview_strains1')
		self.treeview_strains2 = xml.get_widget('treeview_strains2')
		self.button_diff_matrix = xml.get_widget('button_diff_matrix')
		self.canvas_diff_matrix = xml.get_widget('canvas_diff_matrix')
		self.combobox_color_scheme = xml.get_widget('combobox_color_scheme')
		self.radiobutton_allele1 = xml.get_widget('radiobutton_allele1')
		self.radiobutton_allele2 = xml.get_widget('radiobutton_allele2')
		self.radiobutton_het = xml.get_widget('radiobutton_het')
		self.radiobutton_na = xml.get_widget('radiobutton_na')
		self.treeview_diff_details = xml.get_widget('treeview_diff_details')
		self.button_draw_cluster_plot = xml.get_widget('button_draw_cluster_plot')
		
		self.vbox_matplotlib = xml.get_widget('vbox_matplotlib')
		
		# matplotlib canvas
		fig = Figure(figsize=(8,8))
		self.canvas_matplotlib = FigureCanvas(fig)  # a gtk.DrawingArea
		#self.canvas.mpl_connect('button_press_event', self.on_click)
		self.vbox_matplotlib.pack_start(self.canvas_matplotlib)
		
		# matplotlib toolbar
		toolbar = NavigationToolbar(self.canvas_matplotlib, self.app1)
		self.vbox_matplotlib.pack_start(toolbar, False, False)
		
		self.app1_appbar1 = xml.get_widget('app1_appbar1')
		self.app1_appbar1.push('Status Message.')
		self.app1.show_all()
		
		self.app_input = xml.get_widget('app_input')
		self.app_input.connect("destroy", gtk.main_quit)
		
		self.app_input.show_all()
		self.app_vbox_input = xml.get_widget('app_vbox_input')
		self.combobox_qc_class_choice = xml.get_widget('combobox_qc_class_choice')
		self.button_check_class_doc = xml.get_widget('button_check_class_doc')
		self.entry_hostname = xml.get_widget('entry_hostname')
		self.entry_dbname = xml.get_widget('entry_dbname')
		self.entry_schema = xml.get_widget('entry_schema')
		self.entry_1st_input_fname = xml.get_widget('entry_1st_input_fname')
		self.entry_2nd_input_fname = xml.get_widget('entry_2nd_input_fname')
		self.entry_ecotype_table = xml.get_widget('entry_ecotype_table')
		self.entry_ecotype2accession_table = xml.get_widget('entry_ecotype2accession_table')
		self.entry_diff_details_table = xml.get_widget('entry_diff_details_table')
		self.entry_qc_cross_match_table = xml.get_widget('entry_qc_cross_match_table')
		self.entry_latex_output_fname = xml.get_widget('entry_latex_output_fname')
		self.textview_class_doc = xml.get_widget('textview_class_doc')
		#2008-01-21 setting the status in a standalone app bar doesn't work 
		#self.appbar_input = gnome.ui.AppBar(True)
		#self.appbar_input.set_default("Status")
		#self.app_vbox_input.pack_start(self.appbar_input, False, False)
		self.app_input_appbar = xml.get_widget('app_input_appbar')
		self.app_input_appbar.push('Status Message.')
		self.button_input_run = xml.get_widget('button_input_run')
		
		self.qc_class_dict = {0:Cmp250kVs2010.Cmp250kVs2010,
			1:Cmp250kVs149SNP.Cmp250kVs149SNP,
			2:Cmp149SNPVs2010.Cmp149SNPVs2010,
			3:CmpTina2010VsMy2010In250kSNPs.CmpTina2010VsMy2010In250kSNPs}
		self.qc_parent_class_dict = {0:Cmp250kVs2010,
			1:Cmp250kVs149SNP,
			2:Cmp149SNPVs2010,
			3:CmpTina2010VsMy2010In250kSNPs}
		
		self.qc_class_doc = ''
		self.qc_class_ins = None

	def subwindow_hide(self, widget, event, data=None):
		widget.hide()
		return True
	
	def on_combobox_qc_class_choice_changed(self, widget, event=None, data=None):
		class_chosen = self.combobox_qc_class_choice.get_active()
		self.hostname = self.entry_hostname.get_text()
		self.dbname = self.entry_dbname.get_text()
		self.schema = self.entry_schema.get_text()
		import MySQLdb
		self.conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		self.curs = self.conn.cursor()
		self.input_fname1 = self.entry_1st_input_fname.get_text()
		self.input_fname2 = self.entry_2nd_input_fname.get_text()
		self.ecotype_table = self.entry_ecotype_table.get_text()
		self.ecotype2accession_table = self.entry_ecotype2accession_table.get_text()
		self.diff_details_table = self.entry_diff_details_table.get_text()
		self.qc_cross_match_table = self.entry_qc_cross_match_table.get_text()
		self.latex_output_fname = self.entry_latex_output_fname.get_text()
		self.qc_class_doc = self.qc_parent_class_dict[class_chosen].__doc__
		if class_chosen == 2:
			self.qc_class_ins = self.qc_class_dict[class_chosen](self.curs, self.input_fname1, self.input_fname2, self.latex_output_fname, self.ecotype2accession_table, self.diff_details_table, self.qc_cross_match_table)
		else:
			self.qc_class_ins = self.qc_class_dict[class_chosen](self.curs, self.input_fname1, self.input_fname2)

	def on_button_check_class_doc_clicked(self, widget, event=None, data=None):
		self.textbuffer_class_doc = self.textview_class_doc.get_buffer()
		self.textbuffer_class_doc.set_text('')
		if self.qc_class_ins:
			self.textbuffer_class_doc.set_text(self.qc_class_doc)
		else:
			self.textbuffer_class_doc.set_text('No Class Chosen!')
	
	def on_button_input_run_clicked(self, widget, event=None, data=None):
		pass

prog = gnome.program_init('QCVisualizeII', '0.1')
#prog.set_property('app-datadir', '/usr/share')
#prog.set_property('default-icon', '/usr/share/pixmaps/apple-green.png')
instance = QCVisualizeII()
gtk.main()
