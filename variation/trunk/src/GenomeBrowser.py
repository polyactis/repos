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
		self.app1.set_default_size(1000, 800)
		
		self.vbox_matplotlib = xml.get_widget('vbox_matplotlib')
		
		# matplotlib canvas
		fig = Figure(figsize=(8,8))
		self.canvas_matplotlib = FigureCanvas(fig)  # a gtk.DrawingArea
		#self.canvas.mpl_connect('button_press_event', self.on_click)
		self.vbox_matplotlib.pack_start(self.canvas_matplotlib)
		
		# matplotlib toolbar
		toolbar = NavigationToolbar(self.canvas_matplotlib, self.app1)
		self.vbox_matplotlib.pack_start(toolbar, False, False)
		self.app1.show_all()
		self.ax = fig.add_subplot(111)
		self.plot(self.ax, self.canvas_matplotlib, '/tmp/simulate.pvalue')
	
	def plot(self, ax, canvas, input_fname):
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		x_ls = []
		y_ls = []
		for row in reader:
			chr, pos, pvalue = row
			x_ls.append(int(pos))
			y_ls.append(float(pvalue))
		ax.clear()
		ax.plot(x_ls, y_ls, '.')
		canvas.draw()

prog = gnome.program_init('GenomeBrowser', '0.1')
instance = GenomeBrowser()
gtk.main()
