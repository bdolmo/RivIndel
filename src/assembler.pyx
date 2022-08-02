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
    for read in self.reads:
      i = 0
      for i in range(0, len(read)-self.k):
        km1_l = read[i:i+self.k]
        km1_r = read[i+1:i+self.k+1]
        node_l = None
        node_r = None
        if km1_l in nodes:
          node_l = nodes[km1_l]
        else:
          node_l = Node(km1_l)
          nodes[km1_l] = node_l
          G[node_l] = list()
        if km1_r in nodes:
          node_r = nodes[km1_r]
        else:
          node_r = Node(km1_r)
          nodes[km1_r] = node_r
          G[node_r] = list()
        node_r.nin+=1
        node_l.nout+=1
        G[node_l].append(node_r)
    return G

  def eulerian_walk(self) -> str:
    '''
    '''
    start = list(self.g.keys())[0]
    for node in self.g:
      if node.nin < start.nin:
        start = node
    current = start
    contig = current.km1mer

    while len(self.g[current]) > 0:
      next = self.g[current][0]
      del self.g[current][0]
      contig += next.km1mer[-1]
      current = next
    return contig
