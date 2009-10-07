try:
	import numpy as num
	from numpy import *
except:
	try:
		import numarray as num
		from numarray import *
	except:
		import Numeric as num
		from Numeric import *


class CircularQueue:
	'''
	This class implements a data structure for multiple, nested circular queues
	using a numpy array.  The class has several extra useful functions used by NPUTE
	to implement sliding window(s).
	'''
	def __init__(self,Ls,width):
		'''
		Constructor initializes datastructure.
		'''
		height = max(Ls)*2+2
		self.height = height
		self.Ls = Ls
		self.queue = zeros((height,width),uint16)
		self.half = max(Ls)
		self.mid = -1
		
			
	def getEnds(self,i):
		'''
		Returns the bounding rows of the nested queue specified by i.
		Used to calculate mismatch vector over window i.
		'''
		top = self.queue[(self.mid+self.Ls[i])%self.height]
		bottom = self.queue[(self.mid-self.Ls[i]-1)%self.height]
		return top,bottom
		
	def enqueue(self,e):
		'''
		Adds element e to the end of the largest queue.
		'''
		self.incrementMid()
		nextIn = (self.mid+self.half)%self.height
		self.queue[nextIn] = e
		
	def incrementMid(self):
		'''
		Increments the index to the center element of all queues.
		'''
		self.mid = (self.mid+1) % self.height

	def getMid(self):
		'''
		Returns the center element of all queues. Used to remove bias
		when imputing called values during window testing.
		'''
		return self.queue[self.mid]