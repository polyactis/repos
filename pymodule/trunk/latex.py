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