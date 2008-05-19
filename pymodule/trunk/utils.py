import os, sys

def dict_map(dict, ls, type=1):
	"""
	2008-04-03 copied from annot.codense.common
	10-13-05
		add type 2 to return item itself if mapping is not available
	2006-12-21
		add type 3 to extract a smaller map
	2007-05-14
		bug, "if value" could miss 0
	"""
	if type==3:
		new_list = {}	#it's a dictionary
		for item in ls:
			value = dict.get(item)
			if value is not None:
				new_list[item] = value
	else:
		new_list = []
		for item in ls:
			value = dict.get(item)
			if value is not None:
				new_list.append(value)
			elif type==2:
				new_list.append(item)
	
	return new_list

class PassingData(object):
	"""
	05/09/08
		a class to hold any data structure
	"""
	def __init__(self, **keywords):
		"""
		2008-5-12
			add keyword handling
		"""
		for argument_key, argument_value in keywords.iteritems():
			setattr(self, argument_key, argument_value)

def importNumericArray():
	"""
	2008-05-18
		give same numpy types (int, int8 ...) to other numeric modules
	2008-05-18
		add "import array as num"
		should put ImportError in except. but whatever
	2008-05-11
		import whatever available array module
	"""
	try:
		import numpy as num
	except:
		numpy_type2other_ls = ['int', 'int8', 'int16', 'int32', 'int64', 'int128']
		try:
			import numarray as num
		except:
			import Numeric as num
		for numpy_type in self.numpy_type2other_ls:	#make sure it has same type names
			numpy_type_in_other = numpy_type[0].upper() + numpy_type[1:]
			setattr(num, numpy_type, getattr(num, numpy_type_in_other))
	return num

def figureOutDelimiter(input_fname, report=0, delimiter_choice_ls = ['\t', ',']):
	"""
	2008-05-12
		try tab first
	"""
	if report:
		import sys
		sys.stderr.write("Figuring out delimiter for %s ..."%input_fname)
	inf = open(input_fname)
	line = inf.readline()
	del inf
	delimiter_chosen = None
	for delimiter in delimiter_choice_ls:
		delimiter_count = line.count(delimiter)
		if delimiter_count>0:
			delimiter_chosen = delimiter
			break
	if report:
		sys.stderr.write("Done.\n")
	return delimiter_chosen