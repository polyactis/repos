#!/usr/bin/env python
"""
2008-10-07 simple functions related to matplotlib
"""
import os,sys

def assignMatPlotlibHueColorToLs(name_ls, debug=0):
	"""
	2008-10-07
		assign continuous HSL color spectrum to a list (in order).
		
	"""
	if debug:
		sys.stderr.write("Assigning matplotlib hue color to a list ...")
	import ImageColor
	no_of_names = len(name_ls)
	value_step = 1./(no_of_names-1)	#-1 to cover both minimum and maximum in the spectrum
	name2fc = {}
	for i in range(no_of_names):
		name = name_ls[i]
		value = i*value_step
		hue_value = max(min(int(round((1-value)*255)), 255), 0)	#max(min()) makes sure it's 0-255
		fc = ImageColor.getrgb('hsl(%s'%hue_value+',100%,50%)')#"hsl(hue, saturation%, lightness%)" where hue is the colour given as an
		# angle between 0 and 360 (red=0, green=120, blue=240),
		#saturation is a value between 0% and 100% (gray=0%, full color=100%), and lightness is a value between 0% and 100% (black=0%, normal=50%, white=100%).
		fc = [color_value/255. for color_value in fc]	#matplotlib accepts rgb in [0-1] range
		name2fc[name] = fc
	if debug:
		sys.stderr.write("Done.\n")
	return name2fc

def drawName2FCLegend(ax, name_ls, name2fc=None, shape_type=1, no_face_color=False, no_edge_color=False, title=None, font_size=4, alpha=1):
	"""
	2008-10-08
		add option title and linewidth
	2008-10-07
		draw a legend according to name_ls and colors according to name2fc.
		the shape of the legend is symmetric. It's either circle or square.
	"""
	sys.stderr.write("Drawing name2fc legend  ...")
	import matplotlib
	from matplotlib.patches import Polygon, Circle, Ellipse, Wedge
	import numpy
	
	if name2fc is None:
		name2fc = assignMatPlotlibHueColorToLs(name_ls)
	no_of_names = len(name_ls)
	value_step = min(1./(no_of_names), 0.5)	#can't exceed half of plot
	
	radius = value_step/4.
	center_x_pos = 0.25
	center_y_pos = value_step/2.
	xs = [center_x_pos-radius, center_x_pos+radius, center_x_pos+radius, center_x_pos-radius]	#X-axis value for the 4 points of the rectangle starting from lower left corner. 
	ys = numpy.array([center_y_pos-radius, center_y_pos-radius, center_y_pos+radius, center_y_pos+radius])
	for i in range(no_of_names):
		name = name_ls[i]
		fc = name2fc[name]
		facecolor = fc
		edgecolor = fc
		
		if no_face_color:
			facecolor='w'
		if no_edge_color:
			edgecolor='w'
			linewidth = 0
		else:
			linewidth = matplotlib.rcParams['patch.linewidth']
		if shape_type==1:
			patch = Circle((center_x_pos, center_y_pos), radius=radius, linewidth=linewidth, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
			center_y_pos += value_step
		elif shape_type==2:
			patch = Polygon(zip(xs, ys), linewidth=linewidth, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
			ys += value_step	#increase y-axis
		else:
			patch = Circle((center_x_pos, center_y_pos), radius=radius, linewidth=linewidth, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
			center_y_pos += value_step
		ax.add_patch(patch)
		ax.text(0.75, center_y_pos-value_step, name, horizontalalignment ='left', verticalalignment='center', size=font_size)
	if title:
		ax.set_title(title, fontsize=font_size)
	sys.stderr.write("Done.\n")

if __name__ == '__main__':
	#import pdb
	#pdb.set_trace()
	import pylab
	ax = pylab.gca()
	fig = pylab.gcf()
	name_ls = ['1', '2', '3']
	name2fc = assignMatPlotlibHueColorToLs(name_ls)
	drawName2FCLegend(ax, name_ls, name2fc, shape_type=1, no_face_color=True, no_edge_color=False, font_size=10)
	pylab.show()