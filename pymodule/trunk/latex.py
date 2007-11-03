#!/usr/bin/env python
"""
2007-10-15
	modules to deal with latex output
"""

esc_character_dict = {'_':'\\_'}

def escape_characters(latex_line):
	"""
	2007-10-15
		to escape some characters
	"""
	new_latex_line = ''
	for character in latex_line:
		if character in esc_character_dict:
			new_latex_line += esc_character_dict[character]
		else:
			new_latex_line += character
	return new_latex_line

toString = lambda x: '%s'%x

def outputMatrixInLatexTable(data_matrix, caption, table_label, header_ls, footer=None):
	"""
	2007-10-18
		user longtable, rather than tabular
		data_matrix is a two dimensional list
	2007-10-23
		more complicated handling of header_ls
		header_ls could look like [(2,header1), header2, (3,header3)]
		header1 spans 2 columns, header2 1 column, header3 spans 3 columns
	"""
	no_of_rows = len(data_matrix)
	no_of_cols = len(data_matrix[0])
	latex_to_be_returned = '\\begin{center}\n'
	latex_to_be_returned += '\\begin{longtable}{|%s}\n'%('l|'*no_of_cols)
	caption = escape_characters(caption)
	latex_to_be_returned += '\\caption{%s} \\label{%s}\\\\\n'%(caption, table_label)
	latex_to_be_returned += '\\hline\n'
	table_header_in_latex = []
	for header in header_ls:
		if type(header)==list or type(header)==tuple:
			table_header_in_latex.append('\\multicolumn{%s}{|c|}{%s}'%(header[0], header[1]))
		else:
			table_header_in_latex.append(header)
	header_line = '%s\\\\\n'%('&'.join(table_header_in_latex))
	header_line = escape_characters(header_line)
	latex_to_be_returned += header_line
	latex_to_be_returned += '\\hline\n'
	for data_row in data_matrix:	#table entry starts here. one by one
		data_row = map(toString, data_row)
		
		one_table_entry = '&'.join(data_row)
		one_table_entry += '\\\\\n'
		one_table_entry = escape_characters(one_table_entry)
		latex_to_be_returned += one_table_entry
	latex_to_be_returned += '\\hline\n'
	latex_to_be_returned += '\\end{longtable}\n'
	latex_to_be_returned += '\\end{center}\n'
	return latex_to_be_returned

def outputFigureInLatex(fig_fname, caption, fig_label, need_clearpage=0, need_floatpackage=0, fig_width=1):
	"""
	2007-10-19
	2007-10-29
		add need_clearpage, need_floatpackage, fig_width
	"""
	latex_to_be_returned = ''
	if need_clearpage:
		latex_to_be_returned += '\\clearpage\n'
	if need_floatpackage:
		latex_to_be_returned += '\\begin{figure}[H]\n'
	else:
		latex_to_be_returned += '\\begin{figure}\n'
	latex_to_be_returned += '\\includegraphics[width=%s\\textwidth]{%s}\n'%(fig_width,fig_fname)
	caption = escape_characters(caption)
	latex_to_be_returned += '\\caption{%s} \\label{%s}\n'%(caption, fig_label)
	latex_to_be_returned += '\\end{figure}\n'
	return latex_to_be_returned