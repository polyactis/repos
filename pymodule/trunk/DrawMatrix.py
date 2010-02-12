#!/usr/bin/env python
"""
Examples:
	DrawMatrix.py -i /tmp/149CrossMatch_m3.tsv -x /tmp/149CrossMatch_m3.png -s 5
	
	#font size=10, use custom font, matrix value is non-negative. only some special negative value is not ignored.
	#value -2 is colored as black
	DrawMatrix.py -i ./149CrossMatch_m4_a0.3.tsv -x ./banyan_fs/tmp/149CrossMatch_m4_a0.3.png -s 10 -e ./FreeSerif.ttf -m -u black
	
	#partition rows into 2500-row blocks. Color legend is cut into 10 bands.
	DrawMatrix.py -i ./149CrossMatch_m4_a0.1.tsv -x ./banyan_fs/tmp/149CrossMatch_m4_a0.1.png -s 10 -e ./FreeSerif.ttf -m -u black -o 2500 -t 10
	
	#draw the whole matrix into one figure with small font (size=2) and do NOT draw the cell-separating grid.
	DrawMatrix.py -i ./149CrossMatch_m4_a0.1.tsv -x ./banyan_fs/tmp/149CrossMatch_m4_a0.1.png -s 2 -e ./FreeSerif.ttf -m -u black -o 7000 -c 7000 -t 20 -n
	
Description:
	2007-10-23
	module to draw matrix into an image
	some functions copied form annot.bin.codense.common and variation.src.common
	
	2008-08-29
	DrawMatrix is a class to draw matrix values (>=0) in HSL color. values <0 are deemed as special values.
"""

"""
10-31-05 basic functions to draw images
"""
import sys,os
def get_char_dimension():
	import Image, ImageDraw
	im = Image.new('RGB', (50,50))
	draw = ImageDraw.Draw(im)
	char_dimension = draw.textsize('a')
	del im, draw
	return char_dimension


def get_text_region(text, dimension, rotate=1, foreground=(0,0,255), background=(255,255,255), font=None):
	"""
	10-31-05 add background and foreground
	2006-12-13 add font
	"""
	import Image, ImageDraw
	text_im = Image.new('RGB', dimension, background)
	text_draw = ImageDraw.Draw(text_im)
	if font:
		text_draw.text((0,0), text, font=font, fill=foreground)
	else:
		text_draw.text((0,0), text, fill=foreground)
	box = (0,0,dimension[0], dimension[1])
	text_reg = text_im.crop(box)
	if rotate:
		text_reg = text_reg.transpose(Image.ROTATE_90)	#90 is anti-clockwise
	return text_reg

"""
11-28-05
	draw grid to an Image object.
2007-11-02
	copied from annot.bin.codense.common
"""
def get_font(font_path = '/usr/share/fonts/truetype/freefont/FreeSerif.ttf', font_size=20):
	"""
	2008-08-01
		add font_size option
	"""
	import ImageFont
	font = ImageFont.truetype(font_path, font_size)
	return font

def draw_grid(image_object, draw_object, region_to_draw, x_gap, y_gap, color='black'):
	start_x, start_y, stop_x, stop_y = region_to_draw
	#draw horizontal 1st
	for i in range(start_y, stop_y, y_gap):
		draw_object.line((start_x, i+y_gap, stop_x, i+y_gap), fill=color)
	#draw vertical 2nd
	for i in range(start_x, stop_x, x_gap):
		draw_object.line((i+x_gap, start_y, i+x_gap, stop_y), fill=color)

def drawLegend(matrix_value2label, matrix_value2color, font=None):
	"""
	2008-08-21
		matrix_value2color could be a dictionary or function. if it's a dictionary, turn it into lambda function.
	2007-10-23
	2007-11-02
		add line to get default font
	"""
	sys.stderr.write("Drawing legend ...")
	import Image, ImageDraw
	if type(matrix_value2color)==dict:
		matrix_value2color_func = lambda x: matrix_value2color[x]
	else:
		matrix_value2color_func = matrix_value2color
	if not font:
		font = get_font()
	char_dimension = font.getsize('W')	#W is the the biggest(widest)
	char_width, char_height = char_dimension
	label_ls = []
	matrix_value_ls = []
	for matrix_value, color in matrix_value2color.iteritems():
		matrix_value_ls.append(matrix_value)
		label_ls.append(matrix_value2label[matrix_value])
	max_label_len = max(map(len, label_ls))
	label_dimension = (char_width*max_label_len, char_height)
	x_offset0 = 0
	x_offset1 = char_height	#sample color starts here
	x_offset2 = x_offset1 + 2*char_height	#label starts here
	x_offset3 = x_offset2 + label_dimension[0]	#the margin to the right starts here
	y_offset0 = 0
	y_offset1 = y_offset0 + char_height	#sample color starts here
	y_offset2 = y_offset1 + len(label_ls)*2*char_height-char_height	# len(label_ls)-1 gaps among char_height's
	whole_dimension = (x_offset3+char_height, \
			y_offset2+char_height)
	im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
	draw = ImageDraw.Draw(im)
	matrix_value_ls.sort()
	for i in range(len(matrix_value_ls)):
		matrix_value = matrix_value_ls[i]
		y_offset_upper = y_offset1 + i*2*char_height
		y_offset_lower = y_offset1 + i*2*char_height + char_height
		#draw a sample color for this label
		draw.rectangle((x_offset1, y_offset_upper, x_offset1+char_height, y_offset_lower), fill=matrix_value2color_func(matrix_value))
		
		#draw the label
		label = matrix_value2label[matrix_value]
		text_region = get_text_region(label, label_dimension, rotate=0, font=font)	#no rotate
		box = (x_offset2, y_offset_upper, x_offset3, y_offset_lower)
		im.paste(text_region, box)
	#im = im.rotate(270)
	sys.stderr.write("Done.\n")
	return im


def drawContinousLegend(min_value, max_value, no_of_ticks, value2color, font=None, no_of_bands_per_char_height=5):
	"""
	2008-09-30
		fix a bug when min_value, max_value are all integers and band_value_step=0 if no_of_bands>(max_value-min_value).
	2008-08-21
		deal with the case that tick_index goes out of bound
	2008-08-21
		draw legend for continous values
	"""
	sys.stderr.write("Drawing continous legend ...")
	import Image, ImageDraw
	if type(value2color)==dict:
		value2color_func = lambda x: value2color[x]
	else:
		value2color_func = value2color
	if not font:
		font = get_font()
	char_dimension = font.getsize('W')	#W is the the biggest(widest)
	char_width, char_height = char_dimension
	band_height = int(char_height/no_of_bands_per_char_height)
	no_of_bands = 2*(no_of_ticks-1)*no_of_bands_per_char_height	#this is the number of bands to draw
	min_value = float(min_value)	#2008-09-30	convert to float in case all integer cause band_value_step=0 if no_of_bands>(max_value-min_value). python gives integer output if numerator/denominator are all integer.
	max_value = float(max_value)
	band_value_step = (max_value-min_value)/no_of_bands
	band_value_ls = []
	band_value = max_value
	while band_value >= min_value:
		band_value_ls.append(band_value)
		band_value -= band_value_step
	tick_step = (max_value-min_value)/(no_of_ticks-1)
	tick_value_ls = []
	tick_value = max_value
	while tick_value>=min_value:
		tick_value_ls.append(tick_value)
		tick_value -= tick_step
		
	value_label_ls = []
	tick_index = 0
	max_label_len = 0
	for i in range(len(band_value_ls)):
		band_value = band_value_ls[i]
		label = None
		if tick_index<len(tick_value_ls):	#2008-08-21
			tick_value = tick_value_ls[tick_index]
			if abs(band_value-tick_value)<band_value_step:	#if the tick_value and band_value is close enough, bind them together
				label = '%.2f'%tick_value
				if len(label)>max_label_len:
					max_label_len = len(label)
				tick_index += 1			
		value_label_ls.append((band_value, label))
	
	label_dimension = (char_width*max_label_len, char_height)
	x_offset0 = 0
	x_offset1 = char_height	#sample color starts here
	x_offset2 = x_offset1 + 2*char_height	#label starts here
	x_offset3 = x_offset2 + label_dimension[0]	#the margin to the right starts here
	y_offset0 = 0
	y_offset1 = y_offset0 + char_height	#sample color starts here
	y_offset2 = y_offset1 + len(value_label_ls)*band_height	# len(label_ls)-1 gaps among char_height's
	whole_dimension = (x_offset3+char_height, \
			y_offset2+char_height)
	im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
	draw = ImageDraw.Draw(im)
	for i in range(len(value_label_ls)):
		band_value, label = value_label_ls[i]
		
		y_offset_upper = y_offset1 + i*band_height
		y_offset_lower = y_offset1 + (i+1)*band_height
		#draw a sample color for this label
		draw.rectangle((x_offset1, y_offset_upper, x_offset1+char_height, y_offset_lower), fill=value2color_func(band_value))
		
		if label!=None:
			#draw a line here
			draw.line((x_offset1, y_offset_upper, x_offset1+char_height, y_offset_upper), fill='black')
			#draw the label
			text_region = get_text_region(label, label_dimension, rotate=0, font=font)	#no rotate
			box = (x_offset2, y_offset_upper, x_offset3, y_offset_upper+label_dimension[1])
			im.paste(text_region, box)
	sys.stderr.write("Done.\n")
	return im

def drawMatrix(matrix, matrix_value2color, left_label_ls=[], top_label_ls=[], right_label_ls=[], bottom_label_ls=[], with_grid=0, font=None):
	"""
	2008-08-21
		shift the left_label to the right to stick next to the matrix
	2008-08-21
		matrix_value2color could be a dictionary or function. if it's a dictionary, turn it into lambda function.
		use real font-rendition length as maximum label length
	2007-10-23
		matrix is either Numeric, numarray or numpy array.
		use PIL to draw a matrix
	2007-11-02
		add line to get default font
	"""
	sys.stderr.write("Drawing matrix ...")
	import Image, ImageDraw
	if not font:
		font = get_font()
	if type(matrix_value2color)==dict:
		matrix_value2color_func = lambda x: matrix_value2color[x]
	else:
		matrix_value2color_func = matrix_value2color
	char_dimension = font.getsize('W')	#W is the the biggest(widest)
	char_width, char_height = char_dimension
	
	word_font_length = lambda word: font.getsize(word)[0]
	
	if left_label_ls:	#not empty
		max_left_label_length = max(map(word_font_length, left_label_ls))
	else:
		max_left_label_length = 0
	if top_label_ls:
		max_top_label_length = max(map(word_font_length, top_label_ls))
	else:
		max_top_label_length = 0
	if right_label_ls:
		max_right_label_length = max(map(word_font_length, right_label_ls))
	else:
		max_right_label_length = 0
	if bottom_label_ls:
		max_bottom_label_length = max(map(word_font_length, bottom_label_ls))
	else:
		max_bottom_label_length = 0
	
	left_label_dimension = (max_left_label_length + char_width, char_height)
	top_label_dimension = (max_top_label_length + char_width, char_height)	#need rotation
	right_label_dimension = (max_right_label_length + char_width, char_height)
	bottom_label_dimension = (max_bottom_label_length + char_width, char_height)	#need rotation
	
	x_offset0 = char_width
	x_offset1 = left_label_dimension[0]
	x_offset2 = x_offset1 + matrix.shape[1]*char_height
	y_offset0 = 0
	y_offset1 = top_label_dimension[0]
	y_offset2 = y_offset1 + matrix.shape[0]*char_height
	
	whole_dimension = (x_offset2+right_label_dimension[0], \
			y_offset2+bottom_label_dimension[0])
	im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
	draw = ImageDraw.Draw(im)
	#left label
	if left_label_ls:
		for i in range(len(left_label_ls)):
			left_label = left_label_ls[i]
			_left_label_dimension = font.getsize(left_label)
			text_region = get_text_region(left_label, _left_label_dimension, rotate=0, font=font)	#no rotate
			label_shift = max_left_label_length - _left_label_dimension[0]	#shift the left_label to the right to stick next to the matrix
			box = (x_offset0+label_shift, y_offset1+i*left_label_dimension[1], x_offset1, y_offset1+(i+1)*left_label_dimension[1])
			im.paste(text_region, box)
	
	#draw matrix and top_label_ls and bottom_label_ls
	for i in range(matrix.shape[1]):	#x-axis
		x_offset_left = x_offset1+i*top_label_dimension[1]
		x_offset_right = x_offset1+(i+1)*top_label_dimension[1]
		#draw top_label_ls
		if top_label_ls:
			top_label = top_label_ls[i]
			text_region = get_text_region(top_label, top_label_dimension, rotate=1, font=font)
			box = (x_offset_left, y_offset0, x_offset_right, y_offset1)
			im.paste(text_region, box)
		for j in range(matrix.shape[0]):	#y-axis
			draw.rectangle((x_offset_left, y_offset1+j*left_label_dimension[1], \
					x_offset_right, y_offset1+(j+1)*left_label_dimension[1]), fill=matrix_value2color_func(matrix[j,i]))
		#draw bottom_label_ls
		if bottom_label_ls:
			bottom_label = bottom_label_ls[i]
			text_region = get_text_region(bottom_label, bottom_label_dimension, rotate=1, font=font)
			box = (x_offset_left, y_offset2, x_offset_right, whole_dimension[1])
			im.paste(text_region, box)
	
	#right label
	if right_label_ls:
		for i in range(len(right_label_ls)):
			right_label = right_label_ls[i]
			text_region = get_text_region(right_label, right_label_dimension, rotate=0, font=font)
			box = (x_offset2, y_offset2, x_offset_right, whole_dimension[0])
			im.paste(text_region, box)
	if with_grid:
		draw_grid(im, draw, [x_offset1, y_offset1, x_offset2, y_offset2], char_height, char_height)
	
	#im = im.rotate(270)
	sys.stderr.write("Done.\n")
	return im



"""
2007-03-05
	show an image visualizing SNP data
2007-06-05
	set aspect='auto' in imshow(), the default (pylab.image.rcParams['image.aspect'])='equal', which is bad
2007-11-02 copied from variation.src.common, use pylab
"""
def display_snp_matrix(input_fname, output_fname=None, need_sort=0, need_savefig=0, xlabel='', ylabel=''):
	import csv, Numeric, pylab
	reader = csv.reader(open(input_fname), delimiter='\t')
	header = reader.next()
	data_matrix = []
	for row in reader:
		data_row = row[2:]
		data_row = map(int, data_row)
		data_matrix.append(data_row)
	del reader
	data_matrix.reverse()	#2007-03-06 reverse() due to the imshow()'s y axis starting from bottom
	if need_sort:
		data_matrix.sort()
	data_matrix = Numeric.array(data_matrix)
	
	pylab.clf()
	pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')	#2007-06-05
	pylab.colorbar()
	if xlabel:
		pylab.xticks([data_matrix.shape[1]/2], [xlabel])
	if ylabel:
		pylab.yticks([data_matrix.shape[0]/2], [ylabel])
	if need_savefig:
		pylab.savefig('%s.eps'%output_fname, dpi=300)
		pylab.savefig('%s.svg'%output_fname, dpi=300)
		pylab.savefig('%s.png'%output_fname, dpi=300)
	pylab.show()

def make_snp_matrix_legend(value_ls, label_ls, output_fname=None):
	"""
	2007-10-25
		to pair with display_snp_matrix()
	"""
	import numpy, pylab
	label_ls_copy = label_ls[:]
	data_matrix = numpy.zeros([len(value_ls), 1], numpy.int)
	for i in range(len(value_ls)):
		data_matrix[i,0] = value_ls[i]
	label_ls_copy.reverse()	#pylab put the label on starting from the bottom of the picture
	pylab.clf()
	pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')
	pylab.yticks(range(len(value_ls)), label_ls_copy, fontsize=60, verticalalignment='bottom', horizontalalignment='right')
	pylab.xticks([],[])
	if output_fname:
		pylab.savefig('%s.eps'%output_fname, dpi=300)
		pylab.savefig('%s.svg'%output_fname, dpi=300)
		pylab.savefig('%s.png'%output_fname, dpi=300)
	pylab.show()


def display_matrix_of_component(input_fname, ecotypeid_ls, ecotypeid2pos, output_fname=None, need_sort=0, need_savefig=0):
	"""
	2007-09-20
		display the data from that component
		
	"""
	import csv, Numeric, pylab
	cc_ecotypeid_pos = []
	for ecotypeid in ecotypeid_ls:
		cc_ecotypeid_pos.append(ecotypeid2pos[ecotypeid])
	import numpy
	argsort_index = numpy.argsort(cc_ecotypeid_pos, 0)	#watch it's two dimensional
	ecotypeid2row_index = {}
	ytick_label_ls = []
	cc_size = len(ecotypeid_ls)
	for i in range(cc_size):
		ecotypeid_ls_index = argsort_index[i][1]	#sort based on longitude
		ecotypeid = ecotypeid_ls[ecotypeid_ls_index]
		ecotypeid2row_index[ecotypeid] = i
		ytick_label_ls.append('%s (%.2f, %.2f)'%(ecotypeid, ecotypeid2pos[ecotypeid][0], ecotypeid2pos[ecotypeid][1]))
	reader = csv.reader(open(input_fname), delimiter='\t')
	header = reader.next()
	data_matrix = [0]*cc_size
	for row in reader:
		ecotypeid = int(row[0])
		if ecotypeid in ecotypeid2row_index:
			data_row = row[2:]
			data_row = map(int, data_row)
			data_matrix[ecotypeid2row_index[ecotypeid]] = data_row
	del reader
	data_matrix.reverse()	#2007-03-06 reverse() due to the imshow()'s y axis starting from bottom
	if need_sort:
		data_matrix.sort()
	data_matrix = Numeric.array(data_matrix)
	
	pylab.clf()
	pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')	#2007-06-05
	pylab.colorbar()
	pylab.yticks(range(cc_size), ytick_label_ls)
	if need_savefig:
		pylab.savefig('%s.eps'%output_fname, dpi=300)
		pylab.savefig('%s.svg'%output_fname, dpi=300)
		pylab.savefig('%s.png'%output_fname, dpi=300)
	pylab.show()

def combineTwoImages(im1, im2, font=None):
	"""
	2008-08-21
		combine im1 and im2 horizontally
	"""
	sys.stderr.write("Combining two images ...")
	import Image, ImageDraw
	im1_size = im1.size
	im2_size = im2.size
	if not font:
		font = get_font()
	char_dimension = font.getsize('W')	#W is the the biggest(widest)
	char_width, char_height = char_dimension
	whole_dimension = (im1_size[0]+im2_size[0]+char_width, \
			max(im1_size[1], im2_size[1]))
	im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
	im1_box = (0, 0, im1_size[0], im1_size[1])
	im.paste(im1, im1_box)
	
	im2_box = (im1_size[0]+char_width, 0, whole_dimension[0], im2_size[1])
	im.paste(im2, im2_box)
	sys.stderr.write("Done.\n")
	return im

class Value2Color(object):
	"""
	2008-08-29
		for out of range values, treat them as NA
	
	2008-08-28
		a class handles conversion between numerical value and color
		initial functions copied from OutputTestResultInMatrix.py
		
		3 special values:
		-3: separator values
		-2: super_value is a special value out of (min_value, max_value) range.
		-1: NA_value is a value that means the data is missing.
		
		
	"""
	max_gray_value = 255
	special_value2color = {-1:(255,255,255),\
						-2:'red',\
						-3:'black'}
	#usually -1 is NA_value, -2 is super_value.
	def value2GrayScale(cls, value, min_value=0., max_value=255.):
		"""
		2008-08-21
			color span is (0,0,0) to (255,255,255).
		"""
		if value in cls.special_value2color:
			return cls.special_value2color[value]
		elif value<min_value or value > max_value:	#out of range, treat them as NA
			return cls.special_value2color[-1]
		else:
			Y = (value-min_value)/(max_value-min_value)*(cls.max_gray_value-0)
			R_value = cls.max_gray_value-int(round(Y))	#the smaller the value is, the higher hue_value is.
			#in (R,G,B) mode, the bigger R/G/B is, the darker the color is
			#R_value = int(Y/math.pow(2,8))
			#G_value = int(Y- R_value*math.pow(2,8))
			return (R_value, R_value, R_value)
	
	value2GrayScale = classmethod(value2GrayScale)
	
	max_hue_value = 255	#In Inkscape, the maximum possible hue value, 255, looks almost same as hue=0. cut off before reaching 255.
	#but it's not the case in PIL.
	def value2HSLcolor(cls, value, min_value=0., max_value=255., treat_above_max_as_NA=True):
		"""
		2009-10-18
			change how to deal with out-of-range value. Before, both <min_value and >max_value are treated as NA.
			Now, <min_value is treated as NA. >max_value is treated as max_value by default unless treat_above_max_as_NA is set to True. 
		2008-08-28
			use Hue-Saturation-Lightness (HSL) color to replace the simple gray gradient represented by (R,G,B)
		"""
		if value in cls.special_value2color:
			return cls.special_value2color[value]
		elif value<min_value:
			return cls.special_value2color[-1]
		elif value>max_value:	#out of range, treat them as NA
			if treat_above_max_as_NA:
				return cls.special_value2color[-1]
			else:
				value = max_value
		#else:
		Y = (value-min_value)/(max_value-min_value)*(cls.max_hue_value-0)
		hue_value = cls.max_hue_value-int(round(Y))	#the smaller the value is, the higher hue_value is.
		#in (R,G,B) mode, the bigger R/G/B is, the darker the color is
		#R_value = int(Y/math.pow(2,8))
		#G_value = int(Y- R_value*math.pow(2,8))
		return "hsl(%s"%(hue_value)+",100%,50%)"
	
	value2HSLcolor = classmethod(value2HSLcolor)


def drawMatrixLegend(data_matrix, left_label_ls=[], top_label_ls=[], min_value=None, max_value=None,\
					min_possible_value=0, max_possible_value=None,\
					font_path='/usr/share/fonts/truetype/freefont/FreeSerif.ttf', font_size=20, \
					no_of_ticks=15, with_grid=1, treat_above_max_as_NA=True):
	"""
	2009-10-18
		add argument treat_above_max_as_NA, either as maximum value or NA
	2008-10-26
		a wrapper to draw both legend and data matrix by calling drawContinousLegend() and drawMatrix()
		if min_value is not given, try to figure it out by numpy.min(). If min_possible_value is not None, that's the lower bound.
		ditto for max_value.
	"""
	import numpy
	font = get_font(font_path, font_size=font_size)	#2008-08-01
	if min_value==None:
		min_value = numpy.min(data_matrix)
		if min_possible_value is not None:
			min_value = max(min_possible_value, min_value)
	
	if max_value ==None:
		max_value = numpy.max(data_matrix)
		if max_possible_value is not None:
			max_value = min(max_possible_value, max_value)
	
	value2color_func = lambda x: Value2Color.value2HSLcolor(x, min_value, max_value, treat_above_max_as_NA=treat_above_max_as_NA)
	im_legend = drawContinousLegend(min_value, max_value, no_of_ticks, value2color_func, font)
	#im.save('%s_legend.png'%self.fig_fname_prefix)
	im = drawMatrix(data_matrix, value2color_func, left_label_ls,\
				top_label_ls, with_grid=with_grid, font=font)
	im = combineTwoImages(im, im_legend, font=font)
	return im

class DrawMatrix(object):
	__doc__ = __doc__
	option_default_dict = {('font_path', 1, ):['/usr/share/fonts/truetype/freefont/FreeSerif.ttf', 'e', 1, 'path of the font used to draw labels'],\
							('font_size', 1, int):[20, 's', 1, 'size of font, which determines the size of the whole figure.'],\
							("input_fname", 1, ): [None, 'i', 1, 'Filename that stores data matrix. 1st two columns are labels for rows. Top row is header.'],\
							('min_value_non_negative', 0, ):[0, 'm', 0, 'whether minimum value must be >=0 (minus value has special meaning), force min_value=0 if data_matrix gives negative min_value.'],\
							("fig_fname", 1, ): [None, 'x', 1, 'File name for the figure'],\
							("no_of_ticks", 1, int): [5, 't', 1, 'Number of ticks on the legend'],\
							("split_legend_and_matrix", 0, ): [0, 'p', 0, 'whether to split legend and matrix into 2 different images or not. Default is to combine them.'],\
							('super_value_color', 0,): ["red", 'u', 1, 'color for matrix value -2, like "black", or "red" etc.' ],\
							("col_step_size", 0, int): [200, 'c', 1, 'partition columns into blocks according to this size.'],\
							("row_step_size", 0, int): [3500, 'o', 1, 'partition rows into blocks according to this size.'],\
							("no_grid", 0, ): [0, 'n', 0, 'toggle to remove the grid on top of the whole 2-D structure'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-08-29
		"""
		from ProcessOptions import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def _drawMatrix(self, data_matrix, row_label_ls, col_label_ls, fig_fname, passParam):
		"""
		2008-09-01
			parameters come in wrapped passParam=PassingData()
		2008-08-29
			split out of run() to be reusable
		"""
		im = drawMatrix(data_matrix, passParam.value2color_func, row_label_ls, col_label_ls, with_grid=1-passParam.no_grid, font=passParam.font)
		if passParam.split_legend_and_matrix:
			im.save(fig_fname)
		else:
			im = combineTwoImages(im, passParam.im_legend, font=passParam.font)
			im.save(fig_fname)
		
	def run(self):
		"""
		2008-09-10
			in case chop the whole figure into blocks, swap col_block_index and row_block_index to make row first, column 2nd
		"""
		from SNP import read_data
		from utils import figureOutDelimiter, PassingData
		delimiter = figureOutDelimiter(self.input_fname)
		print delimiter
		header, row_label_ls1, row_label_ls2, data_matrix = read_data(self.input_fname, matrix_data_type=float, delimiter='\t')
		import numpy
		data_matrix = numpy.array(data_matrix)
		min_value = numpy.min(data_matrix)
		if self.min_value_non_negative and min_value < 0:
			min_value = 0
		max_value = numpy.max(data_matrix)
		font = get_font(self.font_path, font_size=self.font_size)
		Value2Color.special_value2color[-2] = self.super_value_color
		value2color_func = lambda x: Value2Color.value2HSLcolor(x, min_value, max_value)
		im_legend = drawContinousLegend(min_value, max_value, self.no_of_ticks, value2color_func, font)
		
		fig_fname_prefix = os.path.splitext(self.fig_fname)[0]
		if self.split_legend_and_matrix:
			im_legend.save('%s_legend.png'%fig_fname_prefix)
		
		no_of_rows, no_of_cols = data_matrix.shape
		passParam = PassingData(value2color_func=value2color_func, im_legend=im_legend, font=font, \
							split_legend_and_matrix=self.split_legend_and_matrix, no_grid=self.no_grid)
		
		if no_of_cols <= self.col_step_size:
			self._drawMatrix(data_matrix, row_label_ls1, header[2:], self.fig_fname, passParam)
		else:	#split into blocks
			no_of_col_blocks = no_of_cols/self.col_step_size+1
			no_of_row_blocks = no_of_rows/self.row_step_size + 1
			for i in range(no_of_col_blocks):
				col_start_index = i*self.col_step_size
				col_end_index = (i+1)*self.col_step_size
				if col_start_index<no_of_cols:
					for j in range(no_of_row_blocks):
						row_start_index = j*self.row_step_size
						row_end_index = (j+1)*self.row_step_size
						if row_start_index<no_of_rows:
							fig_fname = '%s_%s_%s.png'%(fig_fname_prefix, j, i)	#row first, column 2nd
							self._drawMatrix(data_matrix[row_start_index:row_end_index,col_start_index:col_end_index], row_label_ls1[row_start_index:row_end_index], \
											header[2+col_start_index:2+col_end_index], fig_fname, passParam)
			
if __name__ == '__main__':
	from ProcessOptions import ProcessOptions
	main_class = DrawMatrix
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
	matrix_value2label= {0:'NA', 1:'A', 2:'C', 3:'G', 4:'T'}
	matrix_value2color = {0:(0,0,122), 1:(0,0,255), 2:(0,122,122), 3:(122,122,0), 4:(255,0,0)}
	font_path = '/usr/share/fonts/truetype/freefont/FreeSerif.ttf'
	font_size=20
	import ImageFont
	font = ImageFont.truetype(font_path, font_size)
	im = drawLegend(matrix_value2label, matrix_value2color, font)
	im.save('/tmp/legend.png')
	
	import sys, os, math
	bit_number = math.log(sys.maxint)/math.log(2)
	if bit_number>40:       #64bit
		sys.path.insert(0, os.path.expanduser('~/script64'))
	else:
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
	from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	input_fname = './script/variation/data/justin_data.csv'
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	import numpy
	data_matrix = numpy.array(data_matrix)
	
	from variation.src.common import number2nt
	im = drawMatrix(data_matrix, matrix_value2color, strain_acc_list, header[2:], with_grid=1, font=font)
	im.save('/tmp/justin_data.png')
	"""