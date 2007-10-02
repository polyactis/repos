#!/usr/bin/env python
"""
2007-10-01
	module to partition a graph into cliques.
	it's a greedy maximum clique searching algorithm.
	always search for the maximum first, then remove it from the graph and keep searching
"""

from clique import clique

class PartitionGraphIntoCliques:
	def __init__(self, debug=0):
		self.debug = int(debug)
		self.clique_ls = []
		self.clique_ins = clique()
	
	def subPartition(self, graph):
		self.clique_ins.max_clique(graph)
		max_clique_ls = self.clique_ins.max_clique_ls
		self.clique_ls.append(max_clique_ls)
		graph.delete_nodes_from(max_clique_ls)
		if graph.number_of_nodes()==0:	#graph.size() is different, = number_of_edges()
			return
		else:
			self.subPartition(graph)
	
	def partition(self, graph):
		if self.debug:
			import pdb
			pdb.set_trace()
		graph1 = graph.copy()
		self.clique_ls = []
		self.subPartition(graph1)

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
	PartitionGraphIntoCliques_ins = PartitionGraphIntoCliques(debug=1)
	PartitionGraphIntoCliques_ins.partition(G)
	print "cliques:", PartitionGraphIntoCliques_ins.clique_ls
	
