#!/usr/bin/env python
"""
Example of embedding matplotlib in an application and interacting with
a treeview to store data.  Double click on an entry to update plot
data

"""
import __init__	#used to know the path to this file itself
import os, sys, pygtk
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

import yh_gnome, csv, traceback
from SNP import SNPData, read_data
from utils import figureOutDelimiter
from sets import Set

class DataMatrixGuiXYProbe(gtk.Window):
	"""
	2009-3-13
		migrated from QCVisualize.py. now become a standalone program and able to read data from a file and plot ...
		QCVisualize.py inherits from here
	2008-02-05
		embed it into a bigger gnome app, add more buttons, and change the __init__()
	2008-01-01
		class to visualize the results from QualityControl.py
	"""
	def __init__(self, plot_title='', id_is_strain=1, header=None, strain_acc_list=None, category_list=None, data_matrix=None):
		"""
		2008-01-10
			use a paned window to wrap the scrolledwindow and the canvas
			so that the relative size of canvas to the scrolledwindow could be adjusted by the user.
		"""
		prog = gnome.program_init('DataMatrixGuiXYProbe', '0.1')	#this must be called before any initialization for gnome app
		
		program_path = os.path.dirname(__init__.__file__)	#sys.argv[0])
		xml = gtk.glade.XML(os.path.join(program_path, 'DataMatrixGuiXYProbe.glade'))
		xml.signal_autoconnect(self)
		self.app1 = xml.get_widget("app1")
		self.app1.connect("delete_event", gtk.main_quit)
		self.app1.set_default_size(800, 1000)
		self.app1.set_title(plot_title)
		
		self.plot_title = plot_title
		self.id_is_strain = id_is_strain
		self.header = header
		self.strain_acc_list = strain_acc_list
		self.category_list = category_list
		self.data_matrix = data_matrix
		
		self.column_types = None
		self.column_header = None
		self.column_editable_flag_ls = None
		self.list_2d = None
		
		
		self.column_types = None
		self.list_2d = None
		self.column_header = None
		self.editable_flag_ls = None
		
		self.vbox1 = xml.get_widget("vbox1")
		self.treeview_matrix = xml.get_widget("treeview_matrix")
		
		# matplotlib stuff
		fig = Figure(figsize=(8,8))
		self.canvas = FigureCanvas(fig)  # a gtk.DrawingArea
		self._idClick = self.canvas.mpl_connect('button_press_event', self.on_click)
		self.vpaned1 = xml.get_widget("vpaned1")
		self.vpaned1.add2(self.canvas)
		
		#vbox.pack_start(self.canvas, True, True)
		self.ax = fig.add_subplot(111)
		self.treeview_matrix.connect('row-activated', self.plot_row)
		
		toolbar = NavigationToolbar(self.canvas, self.app1)
		self.vbox1.pack_start(toolbar, False, False)
		
		self.checkbutton_label_dot = xml.get_widget('checkbutton_label_dot')
		self.entry_dot_label_column = xml.get_widget('entry_dot_label_column')
		self.entry_x_column = xml.get_widget('entry_x_column')
		self.entry_y_column = xml.get_widget('entry_y_column')
		self.entry_hist_column = xml.get_widget('entry_hist_column')
		self.entry_no_of_bins = xml.get_widget('entry_no_of_bins')	#2009-5-20
		self.entry_plot_title = xml.get_widget('entry_plot_title')
		self.entry_plot_title.set_text(self.plot_title)
		
		self.filechooserdialog_save = xml.get_widget("filechooserdialog_save")
		self.filechooserdialog_save.connect("delete_event", yh_gnome.subwindow_hide)
		
		self.filechooserdialog_open = xml.get_widget("filechooserdialog_open")
		self.filechooserdialog_open.connect("delete_event", yh_gnome.subwindow_hide)
		
		self.app1_appbar1 = xml.get_widget('app1_appbar1')
		self.app1_appbar1.push('Status Message.')	#import gnome.ui has to be executed.
		
		self.treeview_matrix.connect('cursor-changed', self.update_no_of_selected, self.app1_appbar1)
		self.app1.show_all()
		
		#self.add_events(gdk.BUTTON_PRESS_MASK|gdk.KEY_PRESS_MASK|gdk.KEY_RELEASE_MASK)
	
	def on_click(self, event):
		"""
		2009-3-13
			use (x_lim[1]-x_lim[0])/200. as the resolution for a dot to be called identical to a data point.
			similar for the y_data
		2009-3-13
			deal with checkbutton_label_dot, entry_dot_label_column, entry_x_column, entry_y_column
		2008-01-01
			derived from on_click_row() of QualityControl.py
			reaction when user clicked in the plot
		"""
		# get the x and y coords, flip y from top to bottom
		x, y = event.x, event.y
		to_label_dot = self.checkbutton_label_dot.get_active()
		dot_label_column = int(self.entry_dot_label_column.get_text())
		x_column = int(self.entry_x_column.get_text())
		y_column = int(self.entry_y_column.get_text())
		x_lim = self.ax.get_xlim()
		x_grain_size = (x_lim[1]-x_lim[0])/200.
		y_lim = self.ax.get_ylim()
		y_grain_size = (y_lim[1]-y_lim[0])/200.
		if event.button==1:
			if event.inaxes is not None:
				print 'data coords', event.xdata, event.ydata
				for row in self.list_2d:
					x_data = row[x_column]
					y_data = row[y_column]
					if abs(x_data-event.xdata)<x_grain_size and abs(y_data-event.ydata)<y_grain_size:
						info = row[dot_label_column]
						if to_label_dot:
							self.ax.text(event.xdata, event.ydata, info, size=8)
							self.canvas.draw()
						sys.stderr.write("%s: %s, %s: %s, xy=(%s, %s), info: %s.\n"%(self.column_header[0], row[0], self.column_header[1], row[1], x_data, y_data, info))
	
	def plotXY(self, ax, canvas, liststore, plot_title='', chosen_index_ls=[]):
		"""
		2009-3-13
			rename plot_NA_mismatch_rate to plotXY()
		2008-02-05
			chosen_index => chosen_index_ls
		2007-12-14
		"""
		x_column = int(self.entry_x_column.get_text())
		y_column = int(self.entry_y_column.get_text())
		plot_title = self.entry_plot_title.get_text()
		
		min_x = 1
		min_y = 1
		max_x = 0
		max_y = 0
		
		x_ls = []
		y_ls = []
		x_chosen_ls = []
		y_chosen_ls = []
		from sets import Set
		chosen_index_set = Set(chosen_index_ls)
		for i in range(len(liststore)):
			row = liststore[i]
			x = row[x_column]
			y = row[y_column]
			if x<min_x:
				min_x = x
			if x>max_x:
				max_x = x
			if y<min_y:
				min_y = y
			if y>max_y:
				max_y = y
			if i in chosen_index_set:
				x_chosen_ls.append(x)
				y_chosen_ls.append(y)
			else:
				x_ls.append(x)
				y_ls.append(y)
		ax.clear()
		ax.plot(x_ls, y_ls, '.')
		
		"""
		#diagonal line give a rough feeling about the notion, more NA, worse calling
		diagonal_start = min(min_x, min_y)-0.1
		diagonal_end = max(max_x, max_x)+0.1
		ax.plot([diagonal_start, diagonal_end],[diagonal_start, diagonal_end])
		"""
		if x_chosen_ls and y_chosen_ls:	#highlight
			ax.plot(x_chosen_ls, y_chosen_ls, '.', c='r')
		if plot_title:
			ax.set_title(plot_title)
		ax.set_xlabel(self.column_header[x_column])
		ax.set_ylabel(self.column_header[y_column])
		canvas.draw()
	
	def plot_row(self, treeview, path, view_column):
		if self._idClick==None:
			self._idClick = self.canvas.mpl_connect('button_press_event', self.on_click)
		self.plotXY(self.ax, self.canvas, self.liststore, self.plot_title, path)
	
	def setupColumns(self, treeview):
		"""
		2009-3-13
		"""
		if not getattr(self, 'column_header', None):
			sys.stderr.write("Nothing in columns yet.\n")
			return
		self.liststore = gtk.ListStore(*self.column_types)
		#self.add_columns(self.treeview_matrix)
		yh_gnome.create_columns(self.treeview_matrix, self.column_header, self.editable_flag_ls, self.liststore)
		yh_gnome.fill_treeview(self.treeview_matrix, self.liststore, self.list_2d, reorderable=True)
		self.treeselection = self.treeview_matrix.get_selection()
	
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
		self.plotXY(self.ax, self.canvas, self.liststore, self.plot_title, index_ls)
	
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
		2009-5-20
			get the number of bins from entry_no_of_bins 
		2009-3-13
			draw histogram of specific hist_column
		2008-02-06
		"""
		if not getattr(self, 'column_header', None):
			sys.stderr.write("Nothing in columns yet.\n")
			return
		self.ax.clear()
		self.canvas.mpl_disconnect(self._idClick)	#drop the signal handler
		self._idClick = None	#reset the _idClick
		hist_ls = []
		hist_column = int(self.entry_hist_column.get_text())
		for i in range(len(self.liststore)):
			hist_ls.append(self.liststore[i][hist_column])
		self.ax.set_title("Histogram of %s %s"%(self.plot_title, self.column_header[hist_column]))
		no_of_bins = int(self.entry_no_of_bins.get_text())
		self.ax.hist(hist_ls, no_of_bins)
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
	
	def readInDataToPlot(self, input_fname):
		"""
		2009-5-20
			add the column index into the column header for easy picking
		2009-3-13
			wrap the float conversion part into try...except to report what goes wrong
		2009-3-13
		"""
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		self.column_header=reader.next()
		for i in range(len(self.column_header)):
			self.column_header[i] = '%s %s'%(i, self.column_header[i])
		no_of_cols = len(self.column_header)
		self.column_types = [str]*2 + [float]*(no_of_cols-2)
		self.column_editable_flag_ls = [True, True] + [False]*(no_of_cols-2)
		self.list_2d = []		
		for row in reader:
			float_part = row[2:]
			try:
				float_part = map(float, float_part)
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				traceback.print_exc()
			new_row = row[:2]+float_part
			self.list_2d.append(new_row)
		self.setupColumns(self.treeview_matrix)
		self.plotXY(self.ax, self.canvas, self.liststore, self.plot_title)
	
	def readInRawMatrixData(self, input_fname):
		"""
		2009-3-13
		"""
		delimiter = figureOutDelimiter(input_fname)
		self.header, self.strain_acc_list, self.category_list, self.data_matrix = read_data(input_fname, delimiter=delimiter)
		
	def on_imagemenuitem_open_activate(self, widget, data=None):
		"""
		2009-3-13
		"""
		self.filechooserdialog_open.show_all()
	
	def on_button_fileopen_ok_clicked(self, widget, data=None):
		"""
		2009-3-13
		"""
		input_fname = self.filechooserdialog_open.get_filename()
		self.filechooserdialog_open.hide()
		self.readInDataToPlot(input_fname)
	
	def on_entry_plot_title_change(self, widget, data=None):
		"""
		2009-3-13
			upon any change in the entry_plot_title
		"""
		self.plot_title = self.entry_plot_title.get_text()
	
if __name__ == '__main__':
	prog = gnome.program_init('DataMatrixGuiXYProbe', '0.1')
	instance = DataMatrixGuiXYProbe()
	gtk.main()