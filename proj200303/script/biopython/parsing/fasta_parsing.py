#!/usr/bin/python

import string
from Bio.ParserSupport import AbstractConsumer

class TitleSearchConsumer(AbstractConsumer):

	def __init_ (self):
		pass

	def title(self,line):
		location = string.find(line, "NAP")
		print location
		
		if location != -1:

			print line
