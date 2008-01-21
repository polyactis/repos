#!/usr/bin/python2.2

import Bio.Clustalw
import Bio.Align.AlignInfo
from Bio.Alphabet import IUPAC
from sys import *

if len(argv) == 2:
    threshold=40.0
else:
    threshold=argv[2]

align = Bio.Clustalw.parse_file(argv[1], alphabet=IUPAC.protein)
alig_len = align.get_alignment_length()
align_info = Bio.Align.AlignInfo.SummaryInfo(align)
ref_seq = align.get_seq_by_num(0)
pssm = align_info.pos_specific_score_matrix(ref_seq, chars_to_ignore = ['X'])
max = len(align_info.get_column(0))

# ----------------------
y = []
for pos in xrange(alig_len):
    max_percent = 0
    for letter in pssm[pos].keys():
        percent = (pssm[pos][letter] / max) * 100.0
        if letter == 'C' and percent > max_percent:
            max_percent = percent
    y.append(max_percent)


# ---------------------------------------------------
#    plot of Cys positions
#

from Numeric import *
from Tkinter import *
import Pmw

# y should have been build before, and must be a tuple
vector_y=tuple(y) 
vector_x = tuple(arange(len(y)))

root = Tk()
frame = Frame(root)
frame.pack()
g = Pmw.Blt.Graph(frame)
g.pack( expand=1, fill='both' )
g.line_create( "Cys sites", xdata=vector_x, ydata=vector_y )
g.configure(width=1000)
g.configure(height=500)
g.element_configure('Cys sites', symbol='none')
g.axis_configure('x', stepsize=100)

