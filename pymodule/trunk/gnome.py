"""
2008-01-21
	module for some custom gnome functions
"""

def foreach_cb(model, path, iter, pathlist):
	"""
	2008-01-21 copied from annot.bin.codense.common
	04-17-05
		used in gui listview, pathfinding.
	"""
	pathlist.append(path)	
	
def create_columns(treeview, label_list):
	"""
	2008-01-21 copied from annot.bin.codense.common
	04-17-05
		create columns in the treeview in the first refresh
	04-21-05
		remove the old columns and reset the model of treeview
	"""
	import gtk
	tvcolumn_dict = {}
	cell_dict = {}
	#remove old columns
	old_column_list = treeview.get_columns()
	for column in old_column_list:
		treeview.remove_column(column)
	treeview.set_model()
		
	for i in range(len(label_list)):
		tvcolumn_dict[i] = gtk.TreeViewColumn(label_list[i])	# create the TreeViewColumn to display the data
		treeview.append_column(tvcolumn_dict[i])	# add tvcolumn to treeview
		cell_dict[i] = gtk.CellRendererText()	# create a CellRendererText to render the data
		tvcolumn_dict[i].pack_start(cell_dict[i], True)	# add the cell to the tvcolumn and allow it to expand
		# set the cell "text" attribute to column 0 - retrieve text
		# from that column in liststore
		tvcolumn_dict[i].add_attribute(cell_dict[i], 'text', i)
		tvcolumn_dict[i].set_sort_column_id(i)	# Allow sorting on the column

def fill_treeview(treeview, liststore, list_2d, reorderable=True):
	"""
	2008-01-21 copied from annot.bin.codense.common
	04-17-05
	05-20-05
		expand the data in list_2d if it's too short
	"""
	import gtk
	length_of_treeview = len(treeview.get_columns())
	for ls in list_2d:
		data = ls[:]	#copy the list to avoid change the content in ls, 'data=ls' changes the content of ls
		for i in range(length_of_treeview-len(data)):
			data.append('')
		liststore.append(data)
	# set the TreeView mode to be liststore
	treeview.set_model(liststore)

	if reorderable:
		for i in range(len(list_2d[0])):
			# make it searchable
			treeview.set_search_column(i)
		
		# Allow drag and drop reordering of rows
		treeview.set_reorderable(True)
	#setting the selection mode
	treeselection = treeview.get_selection()
	treeselection.set_mode(gtk.SELECTION_MULTIPLE)

class Dummy_File:
	"""
	2008-01-24
		copied from http://www.daa.com.au/pipermail/pygtk/attachments/20031004/dbb34c38/Py_Shell.py
	"""
	def __init__(self, buffer):
		"""Implements a file-like object for redirect the stream to the buffer"""
		self.buffer = buffer
	
	def write(self, text):
		"""Write text into the buffer and apply self.tag"""
		iter=self.buffer.get_end_iter()
		self.buffer.insert(iter,text)
	
	def writelines(self, l):
		map(self.write, l)
	
	def flush(self):
		pass
	
	def isatty(self):
		return 1
