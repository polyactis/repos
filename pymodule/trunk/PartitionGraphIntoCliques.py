#!/usr/bin/env python
"""
2007-10-01
	module to partition a graph into cliques.
	it's a greedy maximum clique searching algorithm.
	always search for the maximum first, then remove it from the graph and keep searching
"""

from clique import clique
from networkx import cliques
import numpy

class PartitionGraphIntoCliques:
	def __init__(self, algorithm_type=1, debug=0):
		"""
		2007-10-02
			add algorithm_type
		"""
		self.algorithm_type = int(algorithm_type)
		self.debug = int(debug)
		
		self.whichSubPartition = {1: self.subPartition,
			2: self.subPartitionNx}
		self.clique_ls = []
		self.clique_ins = clique()
	
	def subPartition(self, graph):
		self.clique_ins.max_clique(graph)
		max_clique = self.clique_ins.max_clique_ls
		self.clique_ls.append(max_clique)
		graph.delete_nodes_from(max_clique)
		if graph.number_of_nodes()==0:	#graph.size() is different, = number_of_edges()
			return
		else:
			self.subPartition(graph)
	
	def subPartitionNx(self, graph):
		max_cover_clique_ls = cliques.find_cliques(graph)
		max_clique_index = numpy.argmax(map(len,max_cover_clique_ls))
		max_clique = max_cover_clique_ls[max_clique_index]
		self.clique_ls.append(max_clique)
		graph.delete_nodes_from(max_clique)
		if graph.number_of_nodes()==0:	#graph.size() is different, = number_of_edges()
			return
		else:
			self.subPartitionNx(graph)
	
	def partition(self, graph):
		if self.debug:
			import pdb
			pdb.set_trace()
		graph1 = graph.copy()	# to avoid getting the original graph emptified.
		self.clique_ls = []
		self.whichSubPartition[self.algorithm_type](graph1)

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
	PartitionGraphIntoCliques_ins = PartitionGraphIntoCliques(debug=0)
	PartitionGraphIntoCliques_ins.partition(G)
	print "cliques:", PartitionGraphIntoCliques_ins.clique_ls
	
	"""
	2007-10-02 a testing case that the clique algorithm from clique.py lasts forever
	import sys, os
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/variation/src')))
	from misc import construct_graph_out_of_identity_pair
	input_fname = os.path.expanduser('~/script/variation/stock20070919/data_d110_c0_5.tsv')
	import cPickle
	inf = open('%s_identity_pair_ls'%os.path.splitext(input_fname)[0], 'r')
	identity_pair_ls = cPickle.load(inf)
	del inf
	strain_iden_g = construct_graph_out_of_identity_pair(identity_pair_ls)
	import networkx as nx
	strain_iden_g_cc = nx.connected_components(strain_iden_g)
	g_cc1_sg = strain_iden_g.subgraph(strain_iden_g_cc[1])
	PartitionGraphIntoCliques_ins.partition(g_cc1_sg)
	print "cliques:", PartitionGraphIntoCliques_ins.clique_ls
	"""