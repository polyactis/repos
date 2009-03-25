#!/usr/bin/env python
"""
Example of embedding matplotlib in an application and interacting with
a treeview to store data.  Double click on an entry to update plot
data

"""
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

import yh_gnome

from sets import Set

from DataMatrixGuiXYProbe import DataMatrixGuiXYProbe

class QCVisualize(DataMatrixGuiXYProbe):
	"""
	2009-3-24
		trunk moved to indepedent DataMatrixGuiXYProbe.py and inherits from it
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
		DataMatrixGuiXYProbe.__init__(self, plot_title=plot_title, id_is_strain=id_is_strain)
		
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
				
		self.column_types = [str]*2 + [float]*2 + [int]*5
		self.column_header = ['id', 'id info', 'NA rate', 'mismatch rate', 'no of NAs', 'no of totals', 'no of mismatches', 'no of non NA pairs', 'index in data matrix']
		self.column_editable_flag_ls = [True, True] + [False]*7
		self.list_2d = self.create_list_2d_for_treeview(self.id2NA_mismatch_rate, self.id2info, self.id2index)
		
		self.setupColumns(self.treeview_matrix)
		self.plotXY(self.ax, self.canvas, self.liststore, self.plot_title)
	
	def create_list_2d_for_treeview(self, id2NA_mismatch_rate, id2info, id2index):
		"""
		2008-09-18
			value id2NA_mismatch_rate is a list of more than 6 entries, here only need 6 of them
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
			row = [repr(id), info] + NA_mismatch_rate[:6] + [index_in_data_matrix]
			list_2d.append(row)
		return list_2d