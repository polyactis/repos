#!/usr/bin/env python
"""
2008-02-02 extend the artists in matplotlib
"""
import math
from matplotlib.numerix import array, arange, sin, cos, pi, Float, sqrt, \
	 matrixmultiply, sqrt, nonzero, equal, asarray, dot, concatenate
from matplotlib.artist import Artist, setp, kwdocd
from matplotlib.cbook import dedent
from matplotlib.patches import Polygon
from matplotlib.collections import PolyCollection, LineCollection

# these are not available for the object inspector until after the
# class is build so we define an initial set here for the init
# function and they will be overridden after object defn
kwdocd['Patch'] = """\
		  alpha: float
		  animated: [True | False]
		  antialiased or aa: [True | False]
		  clip_box: a matplotlib.transform.Bbox instance
		  clip_on: [True | False]
		  edgecolor or ec: any matplotlib color
		  facecolor or fc: any matplotlib color
		  figure: a matplotlib.figure.Figure instance
		  fill: [True | False]
		  hatch: unknown
		  label: any string
		  linewidth or lw: float
		  lod: [True | False]
		  transform: a matplotlib.transform transformation instance
		  visible: [True | False]
		  zorder: any number
		  """

class Gene(Polygon):
	"""2008-02-04 An artist for gene. block-like. based on matplotlib.patches.FancyArrow"""

	def __init__(self, start_ls, end_ls, y=0, width=0.001, x_offset=0, is_arrow=True, length_includes_head=True, \
		head_width=None, head_length=None, shape='full', overhang=0, \
		head_starts_at_zero=False,**kwargs):
		"""2008-02-02 an artist for gene. block-like. based on matplotlib.patches.FancyArrow
		
		Returns a new Arrow.
		
		x_offset: is the value of how far start_ls and end_ls should all be pushed
		
		is_arrow: True if it's an arrow. False if it's just a rectangle
		
		length_includes_head: True if head is counted in calculating the length.

		shape: ['full', 'left', 'right']

		overhang: ratio of the head_length that the arrow is swept back (0 overhang means
		triangular shape).

		head_starts_at_zero: if True, the head starts being drawn at coordinate
		0 instead of ending at coordinate 0.

		Valid kwargs are:
		%(Patch)s

		"""
		if head_width is None:
			head_width = 2 * width
		if head_length is None:
			head_length = 1/(2.0 * abs(end_ls[0]-start_ls[0]))

		distance = abs(end_ls[-1]-start_ls[0])
		if length_includes_head:
			length=distance
		else:
			length=distance+head_length
		
		no_of_blocks = len(start_ls)
		if not distance:
			verts = [] #display nothing if empty
		else:
			"""
			start by drawing horizontal arrow, point (tip of the arrow) at (0,0). the whole arrow sticks on the x-axis (<=0 part)
			Notice: the XY -coordination is not the usual math coordination system. it's canvas coordination. (0,0) is top left. horizontal is y-axis. vertical is x-axis.
			start from the last block, reversely to the 1st block
			"""
			inc_block_length = abs(end_ls[-1]-start_ls[-1])
			hw, hl = head_width, head_length
			if is_arrow:
				left_half_arrow_ls = [
    				[0.0,0.0],				  #tip
    				[-hl, -hw/2.0],			 #leftmost
    				[-hl*(1-overhang), -width/2.0], #meets stem
    				[-inc_block_length, -width/2.0]		 #bottom left
    			]
			else:	#this is just a rectangle
				left_half_arrow_ls = [
    				[0.0,0.0],				  #tip
    				[0.0, -width/2.0],		#left to the tip
    				[-inc_block_length, -width/2.0]		 #bottom left
    			]
			if no_of_blocks==1:	#only one block, seal it
				left_half_arrow_ls.append([-inc_block_length, 0])
			elif no_of_blocks>1:	#more than one block
				left_half_arrow_ls.append([-inc_block_length, -width/6.0])   #leave if open
				for i in range(1, no_of_blocks):
					block_index = -i-1	#backwards
					gap = abs(start_ls[block_index+1]-end_ls[block_index])
					inc_block_start = inc_block_length + gap
					this_block_length = abs(end_ls[block_index]-start_ls[block_index])
					inc_block_length = inc_block_length + gap + this_block_length
					if i!=no_of_blocks-1:	#don't seal it
						left_half_arrow_ls.append([-inc_block_start, -width/6.0])
						left_half_arrow_ls.append([-inc_block_start, -width/2.0])
						left_half_arrow_ls.append([-inc_block_length, -width/2.0])
						left_half_arrow_ls.append([-inc_block_length, -width/6.0])
					else:	#seal it
						left_half_arrow_ls.append([-inc_block_start, -width/6.0])
						left_half_arrow_ls.append([-inc_block_start, -width/2.0])
						left_half_arrow_ls.append([-inc_block_length, -width/2.0])
						left_half_arrow_ls.append([-inc_block_length, 0])
			left_half_arrow = array(left_half_arrow_ls)
			#if we're not including the head, shift up by head length
			if not length_includes_head:
				left_half_arrow += [head_length, 0]
			#if the head starts at 0, shift up by another head length
			if head_starts_at_zero:
				left_half_arrow += [head_length/2.0, 0]
			#figure out the shape, and complete accordingly
			if shape == 'left':
				coords = left_half_arrow
			else:
				right_half_arrow = left_half_arrow*[1,-1]
				if shape == 'right':
					coords = right_half_arrow
				elif shape == 'full':
					coords=concatenate([left_half_arrow,right_half_arrow[::-1]])
				else:
					raise ValueError, "Got unknown shape: %s" % shape
			dx = end_ls[-1]-start_ls[0]
			x = start_ls[0] + x_offset
			dy = 0
			cx = float(dx)/distance
			sx = float(dy)/distance
			M = array([[cx, sx],[-sx,cx]])
			verts = matrixmultiply(coords, M) + (x+dx, y+dy)

		Polygon.__init__(self, map(tuple, verts), **kwargs)
	__init__.__doc__ = dedent(__init__.__doc__) % kwdocd

class ExonCollection(PolyCollection):
	"""
	2008-09-23	test PolyCollection, only draw exons
	"""
	def __init__(self, start_ls, end_ls, y=0, width=0.001, x_offset=0, is_arrow=True, length_includes_head=True, \
		head_width=None, head_length=None, shape='full', overhang=0, \
		head_starts_at_zero=False,**kwargs):
		if head_width is None:
			head_width = 2 * width
		if head_length is None:
			head_length = 1/(2.0 * abs(end_ls[0]-start_ls[0]))

		distance = abs(end_ls[-1]-start_ls[0])
		if length_includes_head:
			length=distance
		else:
			length=distance+head_length
		
		no_of_blocks = len(start_ls)
		if not distance:
			verts_ls = [] #display nothing if empty
		else:
			"""
			start by drawing horizontal arrow, point (tip of the arrow) at (0,0). the whole arrow sticks on the x-axis (<=0 part)
			Notice: the XY -coordination is not the usual math coordination system. it's canvas coordination. (0,0) is top left. horizontal is y-axis. vertical is x-axis.
			start from the last block, reversely to the 1st block
			"""
			verts_ls = []
			for i in range(no_of_blocks):
				block_start = start_ls[i]
				block_end = end_ls[i]
				verts = [[block_start, width/2.0+y], [block_start, -width/2.0+y], [block_end, -width/2.0+y], [block_end, width/2.0+y]]
				
				verts_ls.append(verts)
		PolyCollection.__init__(self, verts_ls, **kwargs)

import matplotlib as mpl
import matplotlib.colors as _colors # avoid conflict with kwarg
import matplotlib.transforms as transforms
import matplotlib.cbook as cbook
class ExonIntronCollection(PolyCollection, LineCollection):
	"""
	2008-10-01
		it now can handle picker_event and alpha on both boxes and lines
	2008-09-23	test PolyCollection, draw both exons and introns
	"""
	def __init__(self, annotated_box_ls, y=0, width=0.001, linewidths = None, box_line_widths = None,
				 colors	= None,
				 antialiaseds = None,
				 linestyle = 'solid',
				 offsets = None,
				 transOffset = None,#transforms.identity_transform(),
				 norm = None,
				 cmap = None,
				 picker=False, pickradius=5, alpha=None, strand=None, facecolors=None, edgecolors=None, **kwargs):
		"""
		2008-10-01
			block_type='exon' has both facecolor and edgecolor red.
			block_type='intron' has no facecolor and black edgecolor.
			arguments edgecolors and facecolors in __init__() could overwrite the ones that generated based on the box_type.
		2008-10-01
			the annotated_box_ls is a list of (start, stop, box_type, is_translated, ...)
			box_type = 'exon' or 'intron' and is_translated= 0 or 1. They are optional. Default is box_type ='intron', is_translated=0.
		2008-09-26
			linewidths controls the width of lines connecting exons or the arrow line, could be a list of floats of just one float.
			box_line_widths controls the width of the boundaries of exon boxes
		"""
		#handle parameters for LineCollection
		if linewidths is None   :
			linewidths   = (mpl.rcParams['lines.linewidth'], )
		if box_line_widths is None   :
			box_line_widths   = (mpl.rcParams['lines.linewidth'], )
		if colors is None	   :
			colors	   = (mpl.rcParams['lines.color'],)
		if antialiaseds is None :
			antialiaseds = (mpl.rcParams['lines.antialiased'], )
		
		self._colors = [_colors.colorConverter.to_rgba(color, alpha) for color in colors]	#2008-10-26 use to_rgba() instead of to_rgba_list() because the latter is not present in matplotlib 0.98.3
		self._aa = self._get_value(antialiaseds)
		self._lw = self._get_value(linewidths)
		self.set_linestyle(linestyle)
		self._uniform_offsets = None
		if offsets is not None:
			offsets = npy.asarray(offsets)
			if len(offsets.shape) == 1:
				offsets = offsets[npy.newaxis,:]  # Make it Nx2.
		if transOffset is None:
			if offsets is not None:
				self._uniform_offsets = offsets
				offsets = None
			if hasattr(transforms, 'identity_transform'):
				transOffset = transforms.identity_transform()	#matplotlib 0.91.2
			elif hasattr(transforms, 'IdentityTransform'):
				transOffset = transforms.IdentityTransform()	#matplotlib 0.98.3
		self._offsets = offsets
		self._transOffset = transOffset
		self.pickradius = pickradius
		
		#handle parameters for PolyCollection
		verts_ls = []
		segment_ls = []
		no_of_blocks = len(annotated_box_ls)
		verts_colors = []
		_edgecolors = []
		for i in range(no_of_blocks):
			this_block = annotated_box_ls[i]
			block_start, block_end = this_block[:2]
			block_type = 'intron'
			is_translated = 0
			if len(this_block)>=3:
				block_type = this_block[2]
			if len(this_block)>=4:
				is_translated = int(this_block[3])
			if block_type=='exon':
				verts = [[block_start, width/2.0+y], [block_start, -width/2.0+y], [block_end, -width/2.0+y], [block_end, width/2.0+y]]			
				verts_ls.append(verts)
				if is_translated==1:
					verts_colors.append('r')
					_edgecolors.append('r')	#same edge color
				else:
					verts_colors.append('w')	#no face color
					_edgecolors.append('k')
			else:
				segment_ls.append([(block_start, y), (block_end,y)])
		if not facecolors:
			facecolors = verts_colors
		if not edgecolors:
			edgecolors = _edgecolors
		#add an arrow
		box_start = annotated_box_ls[0][0]
		box_end = annotated_box_ls[-1][1]
		if strand is not None:
			arrow_length = abs(box_start-box_end)/4
			try:
				strand = int(strand)
			except:
				pass
			if strand==1 or strand=='plus':
				arrow_offset = -arrow_length
				arrow_end = box_end
			elif strand==-1 or strand=='minus':
				arrow_offset = arrow_length
				arrow_end = box_start
			else:	#can't tell which strand
				arrow_offset = 0
				arrow_end = box_end
			segment_ls.append([(arrow_end+arrow_offset, y+width), (arrow_end,y+width/2.0)])
		self._segments = segment_ls
		self.set_segments(segment_ls)
		self._verts = verts
		#self.set_alpha(alpha)
		PolyCollection.__init__(self, verts_ls, linewidths=box_line_widths, antialiaseds=antialiaseds, edgecolors=edgecolors, **kwargs)
		self._facecolors  = [_colors.colorConverter.to_rgba(facecolor, alpha=alpha) for facecolor in facecolors]	#in PolyCollection.__init__(), it doesn't handle alpha
		self._picker = picker	#to enable picker_event. PolyCollection's ancestor PatchCollection (collection.py) defines a method 'contains()'.
			#Artist.pick() (in artist.py) calls pickable(), which checks self._picker, first to see whether this artist is pickable. then it calls picker() or contains().
		self.update(kwargs)
	
	def draw(self, renderer):
		if not self.get_visible(): return
		renderer.open_group('polycollection')
		transform = self.get_transform()
		if hasattr(self, 'get_transoffset'):
			transoffset = self.get_transoffset()	#matplotlib 0.91.2
		else:
			transoffset = getattr(self, '_transOffset', None)	#matplotlib 0.98.3

		transform.freeze()
		transoffset.freeze()
		self.update_scalarmappable()
		if cbook.is_string_like(self._edgecolors) and self._edgecolors[:2] == 'No':
			self._linewidths = (0,)
			#self._edgecolors = self._facecolors
		renderer.draw_poly_collection(
			self._verts, transform, self.clipbox,
			self._facecolors, self._edgecolors,
			self._linewidths, self._antialiaseds,
			self._offsets,  transoffset)
		
		self.update_scalarmappable()
		#print 'calling renderer draw line collection'
		offsets = self._offsets
		renderer.draw_line_collection(
			self._segments, transform, self.clipbox,
			self._colors, self._lw, self._ls, self._aa, offsets,
			transoffset)
		
		transform.thaw()
		transoffset.thaw()
		renderer.close_group('polycollection')

def on_canvas_pick(event):
	print dir(event)
	import pdb
	pdb.set_trace()


if __name__ == '__main__':
	import pylab
	a = pylab.gca()
	fig = pylab.gcf()
	fig.canvas.mpl_connect('pick_event', on_canvas_pick)
	start_ls = [1,3,5]
	end_ls = [2,3.5,6]
	pylab.scatter(start_ls, end_ls)
	g = Gene(start_ls, end_ls, y=1,  width=0.1, alpha=0.3, facecolor='r', picker=True)
	a.add_artist(g)
	
	#draw a rectangle
	g1 = Gene(start_ls, end_ls, y=1.5,  width=0.1, x_offset=1, is_arrow=False, alpha=0.3, facecolor='r', picker=True)
	a.add_artist(g1)
	annotated_box_ls = [(1, 1.5, 'exon', 0),(1.5, 2, 'exon', 1), (2,3,'intron'), (3, 3.5, 'exon', 1), (3.5, 5, 'intron'), (5,6,'exon', 1)]
	g2 = ExonIntronCollection(annotated_box_ls, y=0.5, strand=1, width=0.1, facecolors=['y','w', 'w'], alpha=0.3, box_line_widths=0.3, picker=True)
	a.add_artist(g2)
	#draw a opposite strand gene
	start_ls.reverse()
	end_ls.reverse()
	#g2 = ExonCollection(end_ls, start_ls, y=0.5,  width=0.1, x_offset=0.5, head_length=1, facecolors=['r','w', 'w'], overhang=0.3)
	g3 = ExonIntronCollection(annotated_box_ls, strand=-1, y=0,  width=0.1, facecolors=['g','w', 'w'], alpha=0.3, picker=True)
	a.add_artist(g3)
	a.set_xlim(min(start_ls), max(end_ls))
	#a.set_ylim(0,2)
	from matplotlib.patches import Circle, Polygon
	from matplotlib.lines import Line2D
	p1 = Polygon(zip(start_ls, end_ls), alpha=0.3, linewidth=0, facecolor=(0,217/255.,1), picker=True)
	a.add_artist(p1)
	#pylab.gray()	#colorbar needs  "You must first set_array for mappable"
	c1 = Circle((1,2), radius=1, facecolor='w', edgecolor='r')
	a.add_artist(c1)
	#pylab.savefig('abc.png', dpi=300)
	#pylab.savefig('abc.svg', dpi=300)
	pylab.show()