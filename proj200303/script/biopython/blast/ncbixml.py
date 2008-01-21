# Copyright 2002 BioLateral Ltd
# peter.maxwell@biolateral.com.au
# Licence: GPL  http://www.gnu.org/copyleft/gpl.html

"""XML formatted Blast Output Parsing.  BlastOutput = xml_parse(file)"""

__version__ = 1.1

import xml.sax.handler


def safeint(s):
	# only required until int/long unification in approx python2.4
	try:
		v = int(s)
	except ValueError:
		v = long(s)
	return v
			 			

class BlastOutput:
	pass
class Iteration:
	pass	
class Hit:
	pass
class Hsp:
	pass
class Statistics:
	pass
class Parameters:
	pass

class BlastHandler(xml.sax.handler.ContentHandler):
	"""Builds a BlastOutput.  Driven by SAX events from an XML parser"""
	
	attribute_rename = {
		'def': 'definition',  # def is a python keyword
		}
		
	# build attribute_types dictionary for all of the non-integer tags.
	attribute_types = {}
	for tag in ['BlastOutput', 'Parameters', 'Iteration', 'Hit', 'Hsp', 'Statistics']:
		attribute_types[tag] = globals()[tag]
	for (conv, tags) in [
			(str, ['program', 'version', 'reference', 'db', 'query-ID', 'query-def', 'matrix',
				'include', 'filter', 'id', 'def', 'accession', 'qseq', 'hseq', 'midline']),
			(float, ['expect', 'bit-score', 'evalue',
				'density', 'eff-space', 'kappa', 'lambda', 'entropy' ]), # these 5 might be ints
			(list, ['iterations', 'hits', 'hsps']),
			(None, ['param', 'stat'])]:
		for tag in tags:		
			attribute_types[tag] = conv
			
	postprocessed_types = [int, float, str]		
					
	def __init__(self):
		xml.sax.handler.ContentHandler.__init__(self)

	def startDocument(self):
		self.stack = [None]
		
	def endDocument(self):
		assert len(self.stack) == 1 and isinstance(self.stack[0], BlastOutput), self.stack
		
	def startElement(self, tag, attrs):		
		if "_" in tag:
			(parent, tag) = tag.split("_")
		conv = self.attribute_types.get(tag, int)
		if conv is None:
			self.stack.append(None)
			self._expectText = 0
		elif conv in self.postprocessed_types:
			self.stack.append('')
			self._expectText = 1
		else:
			self.stack.append(conv())
			self._expectText = 0
				
	def endElement(self, tag):
		if "_" in tag:
			(parent, tag) = tag.split("_")
		else:
			parent = None
		conv = self.attribute_types.get(tag, int)
		value = self.stack.pop()
		if conv is int:
			value = safeint(value)  
		elif conv in self.postprocessed_types:
			value = conv(value)
		top = self.stack[-1]
		if parent:
			tag = tag.replace('-', '_')
			tag = self.attribute_rename.get(tag, tag)
			assert not tag.startswith('__')
			setattr(top, tag, value)
		elif type(top) is type([]):
			top.append(value)
		elif top is None:
			self.stack[-1] = value	
		else:
			raise RuntimeError, tag		
		self._expectText = 0
		
	def characters(self, text):
		if self._expectText:
			self.stack[-1] = self.stack[-1] + text
		else:
			# There must be a way to tell the xml parser to ignore meaningless 
			# whitespace but I don't know it, so discard it here.
			assert not text.strip(), (text, self._locator.getLineNumber())					
			
def xml_parse(f):
	"""Produce a BlastOutput from an ncbi xml blast output"""
	handler = BlastHandler()
	p = xml.sax.parse(f, handler)
	blast = handler.stack.pop()
	# alias blast.iterations[-1].hits as blast.hits
	final_search = blast.iterations[-1]
	for attr in dir(final_search):
		if not attr.startswith('_'):
			setattr(blast, attr, getattr(final_search, attr))
	return blast
