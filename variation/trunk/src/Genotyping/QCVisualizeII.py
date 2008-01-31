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

from variation.src.QualityControl import QualityControl
from variation.src import Cmp250kVs2010, Cmp250kVs149SNP, Cmp149SNPVs2010, CmpTina2010VsMy2010In250kSNPs

import pymodule.gnome as yh_gnome
#from pymodule.gnome import fill_treeview, create_columns, foreach_cb

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
		self.app1.set_default_size(1400, 800)
		
		self.treeview_strains1 = xml.get_widget('treeview_strains1')
		self.button_diff_matrix = xml.get_widget('button_diff_matrix')
		self.canvas_diff_matrix = xml.get_widget('canvas_diff_matrix')
		canvas_diff_matrix_width = 600
		canvas_diff_matrix_height = 600
		self.canvas_diff_matrix.set_size_request(canvas_diff_matrix_width, canvas_diff_matrix_height)
	        self.canvas_diff_matrix.set_scroll_region(0,0, canvas_diff_matrix_width, canvas_diff_matrix_height)
		self.canvas_diff_matrix.set_data('children_ls',[])	#add a children list (custom data) to canvas_diff_matrix for later cleanup
		self.canvas_diff_matrix.set_data("item_selected", None)	#to designate the item clicked by the user.
		
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
		self.app_input.set_default_size(1000, 800)
		
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
		self.textbuffer_class_doc = self.textview_class_doc.get_buffer()
		self.button_cleanup_output = xml.get_widget('button_cleanup_output')
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
		self.class_chosen = None
		"""
		#2008-01-23 add a callback to monitor IO
		stdout_f = open('/tmp/stdout_f', 'w')
		stderr_f = open('/tmp/stderr_f', 'w')
		sys.stdout = stdout_f
		sys.stderr = stderr_f
		stdout_fr = open('/tmp/stdout_f', 'r')
		stderr_fr = open('/tmp/stderr_f', 'r')
		source_id1 = gobject.io_add_watch(stdout_fr, gobject.IO_IN, self.update_textview_class_doc)
		source_id2 = gobject.io_add_watch(stderr_fr, gobject.IO_IN, self.update_textview_class_doc)
		"""
		#2008-01-24 redirect stdout and stderr to the textbuffer_class_doc
		dummy_file = yh_gnome.Dummy_File(self.textbuffer_class_doc)
		sys.stdout = dummy_file
		sys.stderr = dummy_file
	
	def update_textview_class_doc(self, source, condition):
		"""
		2008-01-24
			dead, doesn't work out. use the Dummy_File instead.
		2008-01-23
			source here is the lower level file descriptor integer and not the Python file object(like sys.stdout)
		"""
		print dir(source)
		source_content = source.read()
		while source_content:
			startiter, enditer = self.textbuffer_class_doc.get_bounds()
			self.textbuffer_class_doc.insert(enditer, source_content)
			#source_content = os.read(source, 20000)
		self.app_input_appbar.push('Output shall be updated.')
	
	def subwindow_hide(self, widget, event, data=None):
		widget.hide()
		return True
	
	def on_combobox_qc_class_choice_changed(self, widget, event=None, data=None):
		self.class_chosen = self.combobox_qc_class_choice.get_active()
		self.qc_class_doc = self.qc_parent_class_dict[self.class_chosen].__doc__
	
	def on_button_check_class_doc_clicked(self, widget, data=None):
		startiter, enditer = self.textbuffer_class_doc.get_bounds()
		#self.textbuffer_class_doc.set_text('')
		if self.class_chosen:
			#self.textbuffer_class_doc.set_text(self.qc_class_doc)
			self.textbuffer_class_doc.insert(enditer, self.qc_class_doc)
			self.app_input_appbar.push('Class doc shown in Output.')
		else:
			self.app_input_appbar.push('Error: No Class Doc cuz No Class Chosen!')
	
	def on_button_cleanup_output_clicked(self, widget, data=None):
		self.textbuffer_class_doc.set_text('')
	
	def on_button_input_run_clicked(self, widget, event=None, data=None):
		if self.class_chosen:
			self.hostname = self.entry_hostname.get_text()
			self.dbname = self.entry_dbname.get_text()
			self.schema = self.entry_schema.get_text()
			import MySQLdb
			self.conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
			self.curs = self.conn.cursor()
			self.input_fname1 = os.path.expanduser(self.entry_1st_input_fname.get_text())
			self.input_fname2 = os.path.expanduser(self.entry_2nd_input_fname.get_text())
			self.ecotype_table = self.entry_ecotype_table.get_text()
			self.ecotype2accession_table = self.entry_ecotype2accession_table.get_text()
			self.diff_details_table = self.entry_diff_details_table.get_text()
			self.qc_cross_match_table = self.entry_qc_cross_match_table.get_text()
			self.latex_output_fname = os.path.expanduser(self.entry_latex_output_fname.get_text())
			if self.class_chosen == 2:
				self.qc_class_ins = self.qc_class_dict[self.class_chosen](self.curs, self.input_fname1, self.input_fname2, self.latex_output_fname, self.ecotype2accession_table, self.diff_details_table, self.qc_cross_match_table)
			else:
				self.qc_class_ins = self.qc_class_dict[self.class_chosen](self.curs, self.input_fname1, self.input_fname2)
			self.qc_class_ins.load_dstruc()
			
			yh_gnome.create_columns(self.treeview_strains1, ['Index', self.qc_class_ins.header1[0], self.qc_class_ins.header2[0], self.qc_class_ins.header2[1]])
			self.strain_index2row_id1 = {}
			self.liststore_strains1 = gtk.ListStore(int, str, str, str)
			strain_acc_list1_2d = [[0, 'ALL', '', '']]
			self.strain_index2row_id1[0] = 'ALL'
			
			row_id22category = {}	#row_id2 's category offers some info on row_id2
			for row_id2, row_index2 in self.qc_class_ins.row_id2row_index2.iteritems():
				row_id22category[row_id2] = self.qc_class_ins.category_list2[row_index2]
			
			strain_index = len(self.strain_index2row_id1) #starting from 1
			for row_id1, row_id2 in self.qc_class_ins.row_id12row_id2.iteritems():
				if row_id2 in row_id22category:
					row_id2_category = row_id22category[row_id2]
				else:
					row_id2_category = ''
				strain_acc_list1_2d.append([strain_index, repr(row_id1), repr(row_id2), row_id2_category])
				self.strain_index2row_id1[strain_index] = row_id1
				strain_index += 1
			
			yh_gnome.fill_treeview(self.treeview_strains1, self.liststore_strains1, strain_acc_list1_2d, reorderable=True)
			
			"""
			#2008-01-24
			yh_gnome.create_columns(self.treeview_strains2, [self.qc_class_ins.header2[0], self.qc_class_ins.header2[1]])
			self.liststore_strains2 = gtk.ListStore(str, str)
			strain_acc_list2_2d = [['ALL', '']]
			for i in range(len(self.qc_class_ins.strain_acc_list2)):
				strain_acc_list2_2d.append([self.qc_class_ins.strain_acc_list2[i], self.qc_class_ins.category_list2[i]])
			yh_gnome.fill_treeview(self.treeview_strains2, self.liststore_strains2, strain_acc_list2_2d, reorderable=True)
			"""
			self.app_input_appbar.push("Check the other window for available strains.")
			self.app1.show_all()
		else:
			self.app_input_appbar.push("Error: Can't Run cuz No Class Chosen!")
	
	def on_button_diff_matrix_clicked(self, widget, data=None):
		"""
		2008-01-21
			figure out the cell-width and cell-height, by adding/deleting a text widget
			set the canvas size, scroll region
			add each widget, set its data, add a callback
		"""
		pathlist_strains1 = []
		treeselection_strains1 = self.treeview_strains1.get_selection()
		treeselection_strains1.selected_foreach(yh_gnome.foreach_cb, pathlist_strains1)
		if len(pathlist_strains1)>0:	#only the first selected one, multiple later
			strain_index = self.liststore_strains1[pathlist_strains1[0][0]][0]
		else:
			strain_index = 0	#default is everything
		
		row_id1 = self.strain_index2row_id1[strain_index]
		if row_id1 == 'ALL':
			row_id1 = -1
		diff_matrix, diff_details_ls, self.diff_code_pair2diff_details_ls = self.qc_class_ins.get_diff_matrix(self.qc_class_ins.data_matrix1, self.qc_class_ins.data_matrix2, self.qc_class_ins.nt_number2diff_matrix_index, self.qc_class_ins.col_id2col_index1, self.qc_class_ins.col_id2col_index2, self.qc_class_ins.col_id12col_id2, self.qc_class_ins.row_id2row_index1, self.qc_class_ins.row_id2row_index2, self.qc_class_ins.row_id12row_id2, row_id=row_id1, need_diff_code_pair_dict=1)
		if diff_matrix!=None:
			for w in self.canvas_diff_matrix.get_data('children_ls'):	#clean up the previous mess
				w.destroy()
				self.canvas_diff_matrix.set_data("item_selected", None)	#reset the item clicked by the user.
			canvas_root = self.canvas_diff_matrix.root()
			no_of_rows, no_of_cols = diff_matrix.shape
			cell_unit_height = self.canvas_diff_matrix.props.height_request/float(no_of_rows)
			cell_unit_width = self.canvas_diff_matrix.props.width_request/float(no_of_cols)
			for i in range(no_of_rows):
				for j in range(no_of_cols):
					w = canvas_root.add(gnomecanvas.CanvasText,text=repr(diff_matrix[i][j]),x=i*cell_unit_height,y=j*cell_unit_width,fill_color='black',anchor=gtk.ANCHOR_W)
					w.set_data('name', (i,j))
					self.canvas_diff_matrix.get_data('children_ls').append(w)
					w.connect("event", self.canvas_diff_matrix_item_event)
		else:
			self.app_input_appbar.push("Error: %s doesn't have corresponding comparison.")
	
	def canvas_diff_matrix_item_event(self, widget, event=None):
		if event.type == gtk.gdk.BUTTON_PRESS:
			if event.button == 1:
				previous_selected_item = self.canvas_diff_matrix.get_data('item_selected')
				if previous_selected_item != None:
					previous_selected_item.set(fill_color='black')	#change the previous selected widget color back to black
				widget.set(fill_color='red')	#change the current selected widget color to red
				self.canvas_diff_matrix.set_data('item_selected', widget)	#remember the current selected widget
				diff_code_pair = widget.get_data('name')
			 	self.app1_appbar1.push(repr(diff_code_pair))
				return True
		elif event.type == gtk.gdk.ENTER_NOTIFY:
			# Make the outline heavy
			widget.set_data('original_weight', widget.get_property('weight'))
			widget.set(weight=1200)
			return True
		elif event.type == gtk.gdk.LEAVE_NOTIFY:
			# Make the outline light
			original_weight = widget.get_data('original_weight')
			widget.set(weight=original_weight)
			return True
		return False

prog = gnome.program_init('QCVisualizeII', '0.1')
#prog.set_property('app-datadir', '/usr/share')
#prog.set_property('default-icon', '/usr/share/pixmaps/apple-green.png')
instance = QCVisualizeII()
gtk.main()
