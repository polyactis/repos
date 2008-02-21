#!/usr/bin/env python
"""
Example of embedding matplotlib in an application and interacting with
a treeview to store data.  Double click on an entry to update plot
data

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

import matplotlib
matplotlib.use('GTKAgg')  # or 'GTK'
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import pymodule.gnome as yh_gnome

from sets import Set

class QCVisualize(gtk.Window):
	"""
	2008-02-05
		embed it into a bigger gnome app, add more buttons, and change the __init__()
	2008-01-01
		class to visualize the results from QualityControl.py
	"""
	def __init__(self, id2NA_mismatch_rate, plot_title='', id2info={}, id2index={}, id_is_strain=1, header=None, strain_acc_list=None, category_list=None, data_matrix=None):
		"""
		2008-01-10
			use a paned window to wrap the scrolledwindow and the canvas
			so that the relative size of canvas to the scrolledwindow could be adjusted by the user.
		"""
		prog = gnome.program_init('QCVisualize', '0.1')	#this must be called before any initialization for gnome app
		
		program_path = os.path.dirname(sys.argv[0])
		xml = gtk.glade.XML(os.path.join(program_path, 'QCVisualize.glade'))
		xml.signal_autoconnect(self)
		self.app1 = xml.get_widget("app1")
		self.app1.connect("delete_event", gtk.main_quit)
		self.app1.set_default_size(800, 1000)
		self.app1.set_title(plot_title)
		
		self.id2NA_mismatch_rate = id2NA_mismatch_rate
		self.plot_title = plot_title
		self.id2info = id2info
		self.id2index = id2index
		self.id_is_strain = id_is_strain
		self.header = header
		self.strain_acc_list = strain_acc_list
		self.category_list = category_list
		self.data_matrix = data_matrix
		
		self.vbox1 = xml.get_widget("vbox1")
		
		types = [str]*2 + [float]*2 + [int]*5
		self.liststore = gtk.ListStore(*types)
		self.list_2d = self.create_list_2d_for_treeview(self.id2NA_mismatch_rate, self.id2info, self.id2index)
		self.treeview_NA_mismatch_rate = xml.get_widget("treeview_NA_mismatch_rate")
		#self.add_columns(self.treeview_NA_mismatch_rate)
		header = ['id', 'id info', 'NA rate', 'mismatch rate', 'no of NAs', 'no of totals', 'no of mismatches', 'no of non NA pairs', 'index in data matrix']
		editable_flag_ls = [True, True] + [False]*7
		yh_gnome.create_columns(self.treeview_NA_mismatch_rate, header, editable_flag_ls, self.liststore)
		yh_gnome.fill_treeview(self.treeview_NA_mismatch_rate, self.liststore, self.list_2d, reorderable=True)
		#self.treeview_NA_mismatch_rate.set_model(self.liststore)
		self.treeselection = self.treeview_NA_mismatch_rate.get_selection()
		#self.treeselection.set_mode(gtk.SELECTION_MULTIPLE)	#already set in fill_treeview()
		
		# matplotlib stuff
		fig = Figure(figsize=(8,8))
		self.canvas = FigureCanvas(fig)  # a gtk.DrawingArea
		self._idClick = self.canvas.mpl_connect('button_press_event', self.on_click)
		self.vpaned1 = xml.get_widget("vpaned1")
		self.vpaned1.add2(self.canvas)
		
		#vbox.pack_start(self.canvas, True, True)
		self.ax = fig.add_subplot(111)
		self.plot_NA_mismatch_rate(self.ax, self.canvas, self.liststore, self.plot_title)
		self.treeview_NA_mismatch_rate.connect('row-activated', self.plot_row)
		
		toolbar = NavigationToolbar(self.canvas, self.app1)
		self.vbox1.pack_start(toolbar, False, False)
		
		self.filechooserdialog_save = xml.get_widget("filechooserdialog_save")
		self.filechooserdialog_save.connect("delete_event", yh_gnome.subwindow_hide)
		
		self.app1_appbar1 = xml.get_widget('app1_appbar1')
		self.app1_appbar1.push('Status Message.')	#import gnome.ui has to be executed.
		
		self.treeview_NA_mismatch_rate.connect('cursor-changed', self.update_no_of_selected, self.app1_appbar1)
		self.app1.show_all()
		
		#self.add_events(gdk.BUTTON_PRESS_MASK|gdk.KEY_PRESS_MASK|gdk.KEY_RELEASE_MASK)
	
	def on_click(self, event):
		"""
		2008-01-01
			derived from on_click_row() of QualityControl.py
			reaction when user clicked in the plot
		"""
		# get the x and y coords, flip y from top to bottom
		x, y = event.x, event.y
		if event.button==1:
			if event.inaxes is not None:
				print 'data coords', event.xdata, event.ydata
				for key, value in self.id2NA_mismatch_rate.iteritems():
					NA_rate, mismatch_rate = value[:2]
					if abs(NA_rate-event.xdata)<0.005 and abs(mismatch_rate-event.ydata)<0.005:
						if key in self.id2info:
							info = self.id2info[key]
						else:
							info = key
						self.ax.text(event.xdata, event.ydata, info, size=8)
						self.canvas.draw()
						print "id: %s, NA_mismatch data: %s, info: %s"%(key, value, info)
	
	def plot_NA_mismatch_rate(self, ax, canvas, liststore, plot_title='', chosen_index_ls=[]):
		"""
		2008-02-05
			chosen_index => chosen_index_ls
		2007-12-14
		"""
		min_NA_rate = 1
		min_mismatch_rate = 1
		max_NA_rate = 0
		max_mismatch_rate = 0
		
		NA_rate_ls = []
		mismatch_rate_ls = []
		NA_rate_chosen_ls = []
		mismatch_rate_chosen_ls = []
		from sets import Set
		chosen_index_set = Set(chosen_index_ls)
		for i in range(len(liststore)):
			row = liststore[i]
			NA_rate = row[2]
			mismatch_rate = row[3]
			if NA_rate<min_NA_rate:
				min_NA_rate = NA_rate
			if NA_rate>max_NA_rate:
				max_NA_rate = NA_rate
			if mismatch_rate<min_mismatch_rate:
				min_mismatch_rate = mismatch_rate
			if mismatch_rate>max_mismatch_rate:
				max_mismatch_rate = mismatch_rate
			if i in chosen_index_set:
				NA_rate_chosen_ls.append(NA_rate)
				mismatch_rate_chosen_ls.append(mismatch_rate)
			else:
				NA_rate_ls.append(NA_rate)
				mismatch_rate_ls.append(mismatch_rate)
		ax.clear()
		ax.plot(NA_rate_ls, mismatch_rate_ls, '.')
		#diagonal line give a rough feeling about the notion, more NA, worse calling
		diagonal_start = min(min_NA_rate, min_mismatch_rate)-0.1
		diagonal_end = max(max_NA_rate, max_mismatch_rate)+0.1
		ax.plot([diagonal_start, diagonal_end],[diagonal_start, diagonal_end])

		if NA_rate_chosen_ls and mismatch_rate_chosen_ls:	#highlight
			ax.plot(NA_rate_chosen_ls, mismatch_rate_chosen_ls, '.', c='r')
		if plot_title:
			ax.set_title(plot_title)
		ax.set_xlabel('NA rate')
		ax.set_ylabel('mismatch rate')
		canvas.draw()
	
	def plot_row(self, treeview, path, view_column):
		if self._idClick==None:
			self._idClick = self.canvas.mpl_connect('button_press_event', self.on_click)
		self.plot_NA_mismatch_rate(self.ax, self.canvas, self.liststore, self.plot_title, path)
	
	def add_columns(self, treeview):
		"""
		2008-02-05
			deprecated, superceded by yh_gnome.create_columns()
		2008-01-01
			add header to the spreadsheet
		"""
		header = ['id', 'id info', 'NA rate', 'mismatch rate', 'no of NAs', 'no of totals', 'no of mismatches', 'no of non NA pairs']
		for i in range(len(header)):
			cellrenderertext = gtk.CellRendererText()
			cellrenderertext.set_property('editable', True)	#set it editable
			column = gtk.TreeViewColumn('%s'%header[i], cellrenderertext, text=i)
			column.set_sort_column_id(i)
			treeview.append_column(column)
			treeview.set_search_column(i)
		treeview.set_reorderable(True)
		
	def create_list_2d_for_treeview(self, id2NA_mismatch_rate, id2info, id2index):
		"""
		2008-02-12
			correct a bug in index_in_data_matrix. it could be 0 and "not index_in_data_matrix" could be come true.
		2008-02-05
			rename from create_model to create_list_2d_for_treeview
			the gtk.ListStore is handled outside.
		2008-01-02
			types = [str]*2 + [float]*2 + [int]*4
		"""
		list_2d = []
		for id, NA_mismatch_rate in id2NA_mismatch_rate.iteritems():
			if id in id2info:
				info = id2info[id]
			else:
				info = ''
			index_in_data_matrix = id2index.get(id)
			if index_in_data_matrix==None:	#2008-02-12 bug
				index_in_data_matrix = -1
			row = [repr(id), info] + NA_mismatch_rate + [index_in_data_matrix]
			list_2d.append(row)
		return list_2d
	
	def on_button_highlight_clicked(self, widget, data=None):
		"""
		2008-02-12
		to update the no_of_selected rows (have to double click a row to change a cursor if it's multiple selection)
		2008-02-05
		"""
		if self._idClick==None:
			self._idClick = self.canvas.mpl_connect('button_press_event', self.on_click)
		pathlist_strains1 = []
		self.treeselection.selected_foreach(yh_gnome.foreach_cb, pathlist_strains1)
		index_ls = []
		for path in pathlist_strains1:
			index_ls.append(path[0])
		self.app1_appbar1.push("%s rows selected."%len(pathlist_strains1))
		self.plot_NA_mismatch_rate(self.ax, self.canvas, self.liststore, self.plot_title, index_ls)
	
	def on_button_save_clicked(self, widget, data=None):
		"""
		2008-02-05
		"""
		self.filechooserdialog_save.show_all()
	
	def on_button_filechooserdialog_cancel_ok_clicked(self, widget, data=None):
		"""
		2008-02-05
		"""
		self.filechooserdialog_save.hide()
	
	def on_button_filechooserdialog_save_ok_clicked(self, widget, data=None):
		"""
		2008-02-12
		to update the no_of_selected rows (have to double click a row to change a cursor if it's multiple selection)
		2008-02-05
		"""
		output_fname = self.filechooserdialog_save.get_filename()
		self.filechooserdialog_save.hide()
		pathlist_strains1 = []
		self.treeselection.selected_foreach(yh_gnome.foreach_cb, pathlist_strains1)
		self.app1_appbar1.push("%s rows selected."%len(pathlist_strains1))
		if self.header and self.strain_acc_list and self.category_list and self.data_matrix:
			selected_index_set = Set()
			for path in pathlist_strains1:
				row = self.liststore[path[0]]
				id = row[0]
				index_in_data_matrix = row[-1]
				selected_index_set.add(index_in_data_matrix)
				if self.id_is_strain:
					id = id[1:-1].split(',')	#id is a tuple of (ecotypeid,duplicate)
					self.strain_acc_list[index_in_data_matrix] = id[0].strip()	#remove extra space
					self.category_list[index_in_data_matrix] = id[1].strip()
				#else:
				#	self.header[index_in_data_matrix+2] = id
			from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
			FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
			if self.id_is_strain:
				rows_to_be_tossed_out = Set(range(len(self.strain_acc_list))) - selected_index_set
				FilterStrainSNPMatrix_instance.write_data_matrix(self.data_matrix, output_fname, self.header, self.strain_acc_list, self.category_list,\
								rows_to_be_tossed_out, cols_to_be_tossed_out=Set(), nt_alphabet=0)
			else:
				cols_to_be_tossed_out = Set(range(len(self.header)-2)) - selected_index_set
				FilterStrainSNPMatrix_instance.write_data_matrix(self.data_matrix, output_fname, self.header, self.strain_acc_list, self.category_list,\
								rows_to_be_tossed_out=Set(), cols_to_be_tossed_out=cols_to_be_tossed_out, nt_alphabet=0)
	
	def show_all(self):
		"""
		2008-02-05
			preserve the old interface. in order not to change anything in plot_col_NA_mismatch_rate() and plot_row_NA_mismatch_rate() of QualityControl.py
		"""
		self.app1.show_all()
	
	def on_button_histogram_clicked(self, widget, data=None):
		"""
		2008-02-06
		"""
		self.ax.clear()
		self.canvas.mpl_disconnect(self._idClick)	#drop the signal handler
		self._idClick = None	#reset the _idClick
		hist_ls = []
		for i in range(len(self.liststore)):
			hist_ls.append(self.liststore[i][3])
		self.ax.set_title("Histogram of %s mismatch rate"%self.plot_title)
		self.ax.hist(hist_ls, 20)
		self.canvas.draw()
	
	def update_no_of_selected(self, treeview, app1_appbar1):
		"""
		2008-02-12
			to update the no_of_selected rows (have to double click a row to change a cursor if it's multiple selection)
		"""
		pathlist_strains1 = []
		self.treeselection.selected_foreach(yh_gnome.foreach_cb, pathlist_strains1)
		app1_appbar1.push("%s rows selected."%len(pathlist_strains1))
		return True