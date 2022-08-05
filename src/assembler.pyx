import os
import sys
import random

class Node():
  '''
  '''
  def __init__(self, km1mer):
    self.km1mer = km1mer
    self.nin = 0
    self.nout = 0

class DeBruijnAssembler():
  '''
  '''
  def __init__(self, reads, k):
    self.reads = reads
    self.k = k
    self.g = self.build_graph()

  def build_graph(self) -> dict:
    '''
      Build a de Bruijn Graph
    '''

    G = dict()
    nodes = dict()

    # unique_reads = set(self.reads)

    for read in self.reads:
      i = 0
      for i in range(0, len(read)-self.k):
        km1_l = read[i:i+self.k]
        km1_r = read[i+1:i+self.k+1]
        node_l = None
        node_r = None
        if km1_l in nodes:
          node_l = nodes[km1_l]
          node_l.nout+=1
        else:
          node_l = Node(km1_l)
          node_l.nout+=1
          nodes[km1_l] = node_l
          G[node_l] = list()
        if km1_r in nodes:
          node_r = nodes[km1_r]
          node_r.nin+=1
        else:
          node_r = Node(km1_r)
          node_r.nin+=1
          nodes[km1_r] = node_r
          G[node_r] = list()
        if node_r not in G[node_l]:
          G[node_l].append(node_r)

    return G

  def eulerian_walk(self):
    '''
    '''
    contig_list = list()

    if len(self.g.keys()) == 0:
      return ""
    start = list(self.g.keys())[0]
    for node in self.g:
      for subnode in self.g[node]:
        if node.nin < start.nin:
          start = node

    current = start
    contig = current.km1mer

    tour = list()
    def _visit(current):
      while len(self.g[current]) > 0:
        dst = self.g[current].pop()
        _visit(dst)
      tour.append(current)

    _visit(current)
    tour = tour[::-1]
    for n in tour:
      if n.nout == 0:
        contig+=n.km1mer[-1]
        contig_list.append(contig)
        contig = ""
      else:
        contig+=n.km1mer[-1]

    return contig_list
