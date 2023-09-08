import numpy as np
import networkx as nx
from common import pmanhatten

class MWPM:

    @staticmethod
    def decode(stabs):
        L = len(stabs)
        ids = np.argwhere(stabs==1) # positions of stabilizers showing "1"

        edges = []
        for i in range(len(ids)-1):
            for j in range(i+1, len(ids)):
                a,b = tuple(ids[i]), tuple(ids[j])
                dist = pmanhatten(a,b,L)
                edges.append( (a,b,dist) )

        G = nx.Graph()
        G.add_weighted_edges_from(edges)
        matchings = nx.algorithms.matching.min_weight_matching(G)
        return matchings
