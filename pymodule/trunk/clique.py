#!/usr/bin/env python
"""
2007-10-01
	module to get the maximum clique
"""
import sys
sys.setrecursionlimit(50000)

class clique:
	def __init__(self, debug=0):
		self.debug = int(debug)
		self.max_clique_size = 0
		self.max_clique_ls = []
	
	def sub_max_clique(self, graph, node_ls, current_clique_size):
		"""
		2007-10-01
			algorithm 1 in Ostergard2002, which is Carraghan1990.
			get some clue about how to save the max clique from Ostergard's cliquer source code
		"""
		no_of_nodes = len(node_ls)
		if no_of_nodes==0:
			if current_clique_size>self.max_clique_size:
				self.max_clique_size = current_clique_size
				self.max_clique_ls = []
				return True
			return False
		return_value = False	#default is False
		while no_of_nodes!=0:
			if current_clique_size+no_of_nodes <= self.max_clique_size:
				return return_value
			v = node_ls.pop(0)	#smallest node
			new_node_ls = []
			for vj in node_ls:	#take the joint set between node_ls and neighbors of v
				if graph.has_edge(v, vj):
					new_node_ls.append(vj)
			if self.sub_max_clique(graph, new_node_ls, current_clique_size+1):
				self.max_clique_ls.append(v)
				if return_value==False:
					return_value = True	#as far as one of the max_clique() returns TRUE, this is successful. in order to trace up all nodes.
			no_of_nodes = len(node_ls)
		return return_value
	
	def order_nodes(self, graph):
		"""
		order the nodes according to degree
		"""
		node_ls = graph.nodes()
		node_degree_ls = graph.degree()
		import numpy
		argsort_node_degree_ls = numpy.argsort(node_degree_ls)
		new_node_ls = []
		for i in argsort_node_degree_ls:
			new_node_ls.append(node_ls[i])
		new_node_ls.reverse()
		return new_node_ls
	
	def init(self):
		self.max_clique_size = 0
		self.max_clique_ls = []
	
	def max_clique(self, graph):
		if self.debug:
			import pdb
			pdb.set_trace()
		self.init()
		node_ls = self.order_nodes(graph)
		self.sub_max_clique(graph, node_ls, 0)

if __name__ == '__main__':
	import networkx as nx
	G=nx.Graph()
	G.add_edge(1,2)
	G.add_edge(2,3)
	G.add_edge(2,4)
	G.add_edge(2,5)
	G.add_edge(2,6)
	G.add_edge(4,6)
	#G.add_edge(7,8)
	#G.add_edge(3,4)
	#G.add_edge(3,6)
	clique_instance = clique(debug=1)
	clique_instance.max_clique(G)
	print 'max_clique_size:', clique_instance.max_clique_size
	print 'max_clique_ls:', clique_instance.max_clique_ls
	print 'nodes:', G.nodes()
	print 'degree:', G.degree()
	nx.draw(G)
	import pylab
	print 'the graph looks like this:'
	pylab.show()
	
	from networkx import cliques
	G_max_clique = cliques.find_cliques(G)
	print 'cliques found  by networkx:', G_max_clique