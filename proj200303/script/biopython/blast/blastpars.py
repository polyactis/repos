#!/usr/bin/python2.2

from Bio.Blast import NCBIStandalone
from Bio.expressions.blast import wublast
from xml.sax import saxutils
parser = wublast.blastn.make_parser(debug_level=0)

outf=open('tmp.xml','w')
parser.setContentHandler(saxutils.XMLGenerator(outf))
parser.parse('tmp2')
wublast.blastn
print wublast.blastn.name


#wublast.main('tmp2')

#bout=open('tmp','r')

#parser=NCBIStandalone.BlastParser()
#record=parser.parse(bout)

#for aln in record.alignments:
#	for hsp in aln.hsps:
#		if hsp.expect<0.00005:
#			print aln.title
#			print aln.length
#			print hsp.query
#			print hsp.sbjct
