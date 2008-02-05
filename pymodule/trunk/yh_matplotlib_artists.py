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



if __name__ == '__main__':
	import pylab
	a = pylab.gca()
	start_ls = [1,3,5]
	end_ls = [2,3.5,6]
	g = Gene(start_ls, end_ls, y=1,  width=0.1, alpha=0.3, facecolor='r', picker=True)
	a.add_artist(g)
	
	#draw a rectangle
	g1 = Gene(start_ls, end_ls, y=1.5,  width=0.1, x_offset=1, is_arrow=False, alpha=0.3, facecolor='r')
	a.add_artist(g1)
	
	#draw a opposite strand gene
	start_ls.reverse()
	end_ls.reverse()
	g2 = Gene(end_ls, start_ls, y=0.5,  width=0.1, x_offset=0.5, head_length=1, alpha=0.3, facecolor='r', overhang=0.3)
	a.add_artist(g2)
	a.set_xlim(min(start_ls), max(end_ls))
	a.set_ylim(0,2)
	
	pylab.show()