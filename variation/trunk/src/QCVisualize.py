#!/usr/bin/env python
"""
Example of embedding matplotlib in an application and interacting with
a treeview to store data.  Double click on an entry to update plot
data

"""
import pygtk
pygtk.require('2.0')
import gtk
from gtk import gdk

import matplotlib
matplotlib.use('GTKAgg')  # or 'GTK'
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar


class QCVisualize(gtk.Window):
	"""
	2008-01-01
		class to visualize the results from QualityControl.py
	"""
	def __init__(self, id2NA_mismatch_rate, plot_title='', id2info={}):
		"""
		2008-01-10
			use a paned window to wrap the scrolledwindow and the canvas
			so that the relative size of canvas to the scrolledwindow could be adjusted by the user.
		"""
		gtk.Window.__init__(self)
		self.id2NA_mismatch_rate = id2NA_mismatch_rate
		self.plot_title = plot_title
		self.id2info = id2info

		self.set_default_size(800, 800)
		self.connect('destroy', lambda win: gtk.main_quit())
		self.set_title('Quality Control Visualize')
		self.set_border_width(4)

		vbox = gtk.VBox(False, 8)
		self.add(vbox)

		label = gtk.Label('Double click a row to highlight it')
		vbox.pack_start(label, False, False)
		
		vpaned = gtk.VPaned()
		vbox.pack_start(vpaned, True, True)
		
		sw = gtk.ScrolledWindow()
		sw.set_shadow_type(gtk.SHADOW_ETCHED_IN)
		sw.set_policy(gtk.POLICY_NEVER,
			          gtk.POLICY_AUTOMATIC)
		vpaned.add1(sw)
		#vbox.pack_start(sw, True, True)

		self.liststore = self.create_model(self.id2NA_mismatch_rate, self.id2info)
		self.treeview = gtk.TreeView(self.liststore)
		self.treeview.set_rules_hint(True)
		self.treeselection = self.treeview.get_selection()
		
		# matplotlib stuff
		fig = Figure(figsize=(8,8))
		self.canvas = FigureCanvas(fig)  # a gtk.DrawingArea
		self.canvas.mpl_connect('button_press_event', self.on_click)
		vpaned.add2(self.canvas)
		
		#vbox.pack_start(self.canvas, True, True)
		self.ax = fig.add_subplot(111)
		
		self.plot_NA_mismatch_rate(self.ax, self.canvas, self.liststore, self.plot_title)
		self.treeview.connect('row-activated', self.plot_row)
		sw.add(self.treeview)

		toolbar = NavigationToolbar(self.canvas, self)
		vbox.pack_start(toolbar, False, False)

		self.add_columns()
		self.add_events(gdk.BUTTON_PRESS_MASK|
				gdk.KEY_PRESS_MASK|
				gdk.KEY_RELEASE_MASK)
	
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
	
	def plot_NA_mismatch_rate(self, ax, canvas, liststore, plot_title='', chosen_index=-1):
		"""
		2007-12-14
		"""
		NA_rate_ls = []
		mismatch_rate_ls = []
		for row in liststore:
			NA_rate = row[2]
			mismatch_rate = row[3]
			NA_rate_ls.append(NA_rate)
			mismatch_rate_ls.append(mismatch_rate)
		ax.clear()
		ax.plot(NA_rate_ls, mismatch_rate_ls, '.')
		#diagonal line give a rough feeling about the notion, more NA, worse calling
		diagonal_start = min(min(NA_rate_ls), min(mismatch_rate_ls))-0.1
		diagonal_end = max(max(NA_rate_ls), max(mismatch_rate_ls))+0.1
		ax.plot([diagonal_start, diagonal_end],[diagonal_start, diagonal_end])

		if chosen_index>=0 or chosen_index<len(liststore):	#highlight
			ax.plot([NA_rate_ls[chosen_index]], [mismatch_rate_ls[chosen_index]], '.', c='r')
		if plot_title:
			ax.set_title(plot_title)
		ax.set_xlabel('NA rate')
		ax.set_ylabel('mismatch rate')
		canvas.draw()
	
	def plot_row(self, treeview, path, view_column):
		ind, = path  # get the index into data
		self.plot_NA_mismatch_rate(self.ax, self.canvas, self.liststore, self.plot_title, ind)
	
	def add_columns(self):
		"""
		2008-01-01
			add header to the spreadsheet
		"""
		header = ['id', 'id info', 'NA rate', 'mismatch rate', 'no of NAs', 'no of totals', 'no of mismatches', 'no of non NA pairs']
		for i in range(len(header)):
			column = gtk.TreeViewColumn('%s'%header[i], gtk.CellRendererText(), text=i)
			column.set_sort_column_id(i)
			self.treeview.append_column(column)
	
	def create_model(self, id2NA_mismatch_rate, id2info):
		"""
		2008-01-02
			types = [str]*2 + [float]*2 + [int]*4
		"""
		types = [str]*2 + [float]*2 + [int]*4
		store = gtk.ListStore(*types)
		for id, NA_mismatch_rate in id2NA_mismatch_rate.iteritems():
			if id in id2info:
				info = id2info[id]
			else:
				info = ''
			row = [repr(id), info] + NA_mismatch_rate
			store.append(row)
		return store
