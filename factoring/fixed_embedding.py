from __future__ import division
from math import sqrt, ceil

import matplotlib.pyplot as plt

import dimod
import dwave_networkx as dnx
import dwave_embedding_utilities as embutil
from dwave.system.samplers import DWaveSampler

sampler = DWaveSampler()
t = 4

nodes_per_cell = t * 2
# edges_per_cell = t * t
m = n = int(ceil(sqrt(ceil(len(sampler.structure.nodelist) / nodes_per_cell))))  # assume square lattice shape
system = dnx.chimera_graph(m, n, t, node_list=sampler.structure.nodelist, edge_list=sampler.structure.edgelist)
c2i = {chimera_index: linear_index for (linear_index, chimera_index) in system.nodes(data='chimera_index')}

bqm = dimod.BinaryQuadraticModel({}, {}, 0, 'SPIN')




bqm.add_variable(c2i[0, 0, 0, 0], 1)
bqm.add_variable(c2i[0, 0, 0, 1], -0.5)
# bqm.add_variable(c2i[0, 0, 0, 2], 0)
# bqm.add_variable(c2i[0, 0, 0, 3], 0)

bqm.add_variable(c2i[0, 0, 1, 0], 0)
# bqm.add_variable(c2i[0, 0, 1, 1], 0)
bqm.add_variable(c2i[0, 0, 1, 2], -0.5)
# bqm.add_variable(c2i[0, 0, 1, 3], 0)


# bqm.add_variable(c2i[0, 1, 0, 0], 0)
# bqm.add_variable(c2i[0, 1, 0, 1], 0)
# bqm.add_variable(c2i[0, 1, 0, 2], 0)
# bqm.add_variable(c2i[0, 1, 0, 3], 0)

# bqm.add_variable(c2i[0, 1, 1, 0], 0)
# bqm.add_variable(c2i[0, 1, 1, 1], 0)
bqm.add_variable(c2i[0, 1, 1, 2], 0)
# bqm.add_variable(c2i[0, 1, 1, 3], 0)


# bqm.add_variable(c2i[0, 2, 0, 0], 0)
# bqm.add_variable(c2i[0, 2, 0, 1], 0)
bqm.add_variable(c2i[0, 2, 0, 2], 1)
bqm.add_variable(c2i[0, 2, 0, 3], -0.5)

bqm.add_variable(c2i[0, 2, 1, 0], 0)
# bqm.add_variable(c2i[0, 2, 1, 1], 0)
bqm.add_variable(c2i[0, 2, 1, 2], -0.5)
# bqm.add_variable(c2i[0, 2, 1, 3], 0)


# bqm.add_variable(c2i[0, 3, 0, 0], 0)
# bqm.add_variable(c2i[0, 3, 0, 1], 0)
# bqm.add_variable(c2i[0, 3, 0, 2], 0)
# bqm.add_variable(c2i[0, 3, 0, 3], 0)

# bqm.add_variable(c2i[0, 3, 1, 0], 0)
# bqm.add_variable(c2i[0, 3, 1, 1], 0)
bqm.add_variable(c2i[0, 3, 1, 2], 0)
# bqm.add_variable(c2i[0, 3, 1, 3], 0)


bqm.add_variable(c2i[0, 4, 0, 0], 1)
bqm.add_variable(c2i[0, 4, 0, 1], -0.5)
# bqm.add_variable(c2i[0, 4, 0, 2], 0)
# bqm.add_variable(c2i[0, 4, 0, 3], 0)

bqm.add_variable(c2i[0, 4, 1, 0], 0)
bqm.add_variable(c2i[0, 4, 1, 1], 0)
bqm.add_variable(c2i[0, 4, 1, 2], -0.5)
# bqm.add_variable(c2i[0, 4, 1, 3], 0)


bqm.add_variable(c2i[1, 0, 0, 0], 0)
bqm.add_variable(c2i[1, 0, 0, 1], 0)
# bqm.add_variable(c2i[1, 0, 0, 2], 0)
# bqm.add_variable(c2i[1, 0, 0, 3], 0)

# bqm.add_variable(c2i[1, 0, 1, 0], 0)
# bqm.add_variable(c2i[1, 0, 1, 1], 0)
bqm.add_variable(c2i[1, 0, 1, 2], 0)
# bqm.add_variable(c2i[1, 0, 1, 3], 0)


bqm.add_variable(c2i[1, 1, 0, 0], 1)
bqm.add_variable(c2i[1, 1, 0, 1], 0)
bqm.add_variable(c2i[1, 1, 0, 2], 1)
# bqm.add_variable(c2i[1, 1, 0, 3], 0)

bqm.add_variable(c2i[1, 1, 1, 0], 0)
bqm.add_variable(c2i[1, 1, 1, 1], 1)
bqm.add_variable(c2i[1, 1, 1, 2], 0)
# bqm.add_variable(c2i[1, 1, 1, 3], 0)


bqm.add_variable(c2i[1, 2, 0, 0], 0)
bqm.add_variable(c2i[1, 2, 0, 1], 0)
bqm.add_variable(c2i[1, 2, 0, 2], 0)
bqm.add_variable(c2i[1, 2, 0, 3], 0)

bqm.add_variable(c2i[1, 2, 1, 0], 0)
bqm.add_variable(c2i[1, 2, 1, 1], 0)
bqm.add_variable(c2i[1, 2, 1, 2], 0)
# bqm.add_variable(c2i[1, 2, 1, 3], 0)


bqm.add_variable(c2i[1, 3, 0, 0], 1)
bqm.add_variable(c2i[1, 3, 0, 1], 0)
bqm.add_variable(c2i[1, 3, 0, 2], 1)
# bqm.add_variable(c2i[1, 3, 0, 3], 0)

bqm.add_variable(c2i[1, 3, 1, 0], 0)
bqm.add_variable(c2i[1, 3, 1, 1], 1)
bqm.add_variable(c2i[1, 3, 1, 2], 0)
# bqm.add_variable(c2i[1, 3, 1, 3], 0)


bqm.add_variable(c2i[1, 4, 0, 0], 0)
bqm.add_variable(c2i[1, 4, 0, 1], 0)
# bqm.add_variable(c2i[1, 4, 0, 2], 0)
# bqm.add_variable(c2i[1, 4, 0, 3], 0)

bqm.add_variable(c2i[1, 4, 1, 0], 0)
# bqm.add_variable(c2i[1, 4, 1, 1], 0)
# bqm.add_variable(c2i[1, 4, 1, 2], 0)
# bqm.add_variable(c2i[1, 4, 1, 3], 0)


bqm.add_variable(c2i[2, 0, 0, 0], 1)
bqm.add_variable(c2i[2, 0, 0, 1], -0.5)
# bqm.add_variable(c2i[2, 0, 0, 2], 0)
# bqm.add_variable(c2i[2, 0, 0, 3], 0)

bqm.add_variable(c2i[2, 0, 1, 0], 0)
bqm.add_variable(c2i[2, 0, 1, 1], -0.5)
# bqm.add_variable(c2i[2, 0, 1, 2], 0)
# bqm.add_variable(c2i[2, 0, 1, 3], 0)


# bqm.add_variable(c2i[2, 1, 0, 0], 0)
# bqm.add_variable(c2i[2, 1, 0, 1], 0)
bqm.add_variable(c2i[2, 1, 0, 2], 0)
# bqm.add_variable(c2i[2, 1, 0, 3], 0)

# bqm.add_variable(c2i[2, 1, 1, 0], 0)
bqm.add_variable(c2i[2, 1, 1, 1], 0)
# bqm.add_variable(c2i[2, 1, 1, 2], 0)
# bqm.add_variable(c2i[2, 1, 1, 3], 0)


bqm.add_variable(c2i[2, 2, 0, 0], 0)
bqm.add_variable(c2i[2, 2, 0, 1], 1)
# bqm.add_variable(c2i[2, 2, 0, 2], 0)
bqm.add_variable(c2i[2, 2, 0, 3], -0.5)

bqm.add_variable(c2i[2, 2, 1, 0], 0)
bqm.add_variable(c2i[2, 2, 1, 1], -0.5)
# bqm.add_variable(c2i[2, 2, 1, 2], 0)
# bqm.add_variable(c2i[2, 2, 1, 3], 0)


# bqm.add_variable(c2i[2, 3, 0, 0], 0)
# bqm.add_variable(c2i[2, 3, 0, 1], 0)
bqm.add_variable(c2i[2, 3, 0, 2], 0)
# bqm.add_variable(c2i[2, 3, 0, 3], 0)

# bqm.add_variable(c2i[2, 3, 1, 0], 0)
bqm.add_variable(c2i[2, 3, 1, 1], 0)
# bqm.add_variable(c2i[2, 3, 1, 2], 0)
# bqm.add_variable(c2i[2, 3, 1, 3], 0)


bqm.add_variable(c2i[2, 4, 0, 0], 1)
bqm.add_variable(c2i[2, 4, 0, 1], -0.5)
# bqm.add_variable(c2i[2, 4, 0, 2], 0)
# bqm.add_variable(c2i[2, 4, 0, 3], 0)

bqm.add_variable(c2i[2, 4, 1, 0], 0)
bqm.add_variable(c2i[2, 4, 1, 1], -0.5)
# bqm.add_variable(c2i[2, 4, 1, 2], 0)
# bqm.add_variable(c2i[2, 4, 1, 3], 0)


bqm.add_variable(c2i[3, 0, 0, 0], 0)
bqm.add_variable(c2i[3, 0, 0, 1], 0)
# bqm.add_variable(c2i[3, 0, 0, 2], 0)
# bqm.add_variable(c2i[3, 0, 0, 3], 0)

# bqm.add_variable(c2i[3, 0, 1, 0], 0)
# bqm.add_variable(c2i[3, 0, 1, 1], 0)
# bqm.add_variable(c2i[3, 0, 1, 2], 0)
bqm.add_variable(c2i[3, 0, 1, 3], 0)


bqm.add_variable(c2i[3, 1, 0, 0], 0)
bqm.add_variable(c2i[3, 1, 0, 1], 0)
bqm.add_variable(c2i[3, 1, 0, 2], 0)
bqm.add_variable(c2i[3, 1, 0, 3], 0)

bqm.add_variable(c2i[3, 1, 1, 0], 0)
bqm.add_variable(c2i[3, 1, 1, 1], 0)
bqm.add_variable(c2i[3, 1, 1, 2], 0)
bqm.add_variable(c2i[3, 1, 1, 3], 0)


bqm.add_variable(c2i[3, 2, 0, 0], 0)
bqm.add_variable(c2i[3, 2, 0, 1], 0)
bqm.add_variable(c2i[3, 2, 0, 2], 0)
bqm.add_variable(c2i[3, 2, 0, 3], 0)

bqm.add_variable(c2i[3, 2, 1, 0], 0)
bqm.add_variable(c2i[3, 2, 1, 1], 0)
# bqm.add_variable(c2i[3, 2, 1, 2], 0)
bqm.add_variable(c2i[3, 2, 1, 3], 0)


bqm.add_variable(c2i[3, 3, 0, 0], 0)
bqm.add_variable(c2i[3, 3, 0, 1], 0)
bqm.add_variable(c2i[3, 3, 0, 2], 0)
bqm.add_variable(c2i[3, 3, 0, 3], 0)

bqm.add_variable(c2i[3, 3, 1, 0], 0)
bqm.add_variable(c2i[3, 3, 1, 1], 0)
bqm.add_variable(c2i[3, 3, 1, 2], 0)
bqm.add_variable(c2i[3, 3, 1, 3], 0)


bqm.add_variable(c2i[3, 4, 0, 0], 0)
bqm.add_variable(c2i[3, 4, 0, 1], 0)
# bqm.add_variable(c2i[3, 4, 0, 2], 0)
# bqm.add_variable(c2i[3, 4, 0, 3], 0)

bqm.add_variable(c2i[3, 4, 1, 0], 0)
# bqm.add_variable(c2i[3, 4, 1, 1], 0)
# bqm.add_variable(c2i[3, 4, 1, 2], 0)
# bqm.add_variable(c2i[3, 4, 1, 3], 0)


bqm.add_variable(c2i[4, 0, 0, 0], 1)
bqm.add_variable(c2i[4, 0, 0, 1], -0.5)
# bqm.add_variable(c2i[4, 0, 0, 2], 0)
# bqm.add_variable(c2i[4, 0, 0, 3], 0)

bqm.add_variable(c2i[4, 0, 1, 0], 0)
bqm.add_variable(c2i[4, 0, 1, 1], -0.5)
# bqm.add_variable(c2i[4, 0, 1, 2], 0)
# bqm.add_variable(c2i[4, 0, 1, 3], 0)


# bqm.add_variable(c2i[4, 1, 0, 0], 0)
# bqm.add_variable(c2i[4, 1, 0, 1], 0)
# bqm.add_variable(c2i[4, 1, 0, 2], 0)
bqm.add_variable(c2i[4, 1, 0, 3], 0)

# bqm.add_variable(c2i[4, 1, 1, 0], 0)
bqm.add_variable(c2i[4, 1, 1, 1], 0)
# bqm.add_variable(c2i[4, 1, 1, 2], 0)
# bqm.add_variable(c2i[4, 1, 1, 3], 0)


bqm.add_variable(c2i[4, 2, 0, 0], 0)
bqm.add_variable(c2i[4, 2, 0, 1], 0)
bqm.add_variable(c2i[4, 2, 0, 2], 1)
bqm.add_variable(c2i[4, 2, 0, 3], -0.5)

bqm.add_variable(c2i[4, 2, 1, 0], 0)
bqm.add_variable(c2i[4, 2, 1, 1], -0.5)
bqm.add_variable(c2i[4, 2, 1, 2], 0)
# bqm.add_variable(c2i[4, 2, 1, 3], 0)


# bqm.add_variable(c2i[4, 3, 0, 0], 0)
# bqm.add_variable(c2i[4, 3, 0, 1], 0)
bqm.add_variable(c2i[4, 3, 0, 2], 0)
bqm.add_variable(c2i[4, 3, 0, 3], 0)

bqm.add_variable(c2i[4, 3, 1, 0], 0)
bqm.add_variable(c2i[4, 3, 1, 1], 0)
# bqm.add_variable(c2i[4, 3, 1, 2], 0)
# bqm.add_variable(c2i[4, 3, 1, 3], 0)


bqm.add_variable(c2i[4, 4, 0, 0], 1)
bqm.add_variable(c2i[4, 4, 0, 1], -0.5)
# bqm.add_variable(c2i[4, 4, 0, 2], 0)
# bqm.add_variable(c2i[4, 4, 0, 3], 0)

bqm.add_variable(c2i[4, 4, 1, 0], 0)
bqm.add_variable(c2i[4, 4, 1, 1], -0.5)
# bqm.add_variable(c2i[4, 4, 1, 2], 0)
# bqm.add_variable(c2i[4, 4, 1, 3], 0)


bqm.add_variable(c2i[5, 0, 0, 0], 0)
# bqm.add_variable(c2i[5, 0, 0, 1], 0)
# bqm.add_variable(c2i[5, 0, 0, 2], 0)
# bqm.add_variable(c2i[5, 0, 0, 3], 0)

# bqm.add_variable(c2i[5, 0, 1, 0], 0)
# bqm.add_variable(c2i[5, 0, 1, 1], 0)
# bqm.add_variable(c2i[5, 0, 1, 2], 0)
bqm.add_variable(c2i[5, 0, 1, 3], 0)


bqm.add_variable(c2i[5, 1, 0, 0], 0)
bqm.add_variable(c2i[5, 1, 0, 1], 0)
bqm.add_variable(c2i[5, 1, 0, 2], 0)
bqm.add_variable(c2i[5, 1, 0, 3], 0)

bqm.add_variable(c2i[5, 1, 1, 0], 0)
bqm.add_variable(c2i[5, 1, 1, 1], 0)
bqm.add_variable(c2i[5, 1, 1, 2], 0)
bqm.add_variable(c2i[5, 1, 1, 3], 0)


bqm.add_variable(c2i[5, 2, 0, 0], 0)
bqm.add_variable(c2i[5, 2, 0, 1], 0)
# bqm.add_variable(c2i[5, 2, 0, 2], 0)
# bqm.add_variable(c2i[5, 2, 0, 3], 0)

bqm.add_variable(c2i[5, 2, 1, 0], 0)
# bqm.add_variable(c2i[5, 2, 1, 1], 0)
bqm.add_variable(c2i[5, 2, 1, 2], 0)
# bqm.add_variable(c2i[5, 2, 1, 3], 0)


bqm.add_variable(c2i[5, 3, 0, 0], 1)
bqm.add_variable(c2i[5, 3, 0, 1], 0)
bqm.add_variable(c2i[5, 3, 0, 2], 1)
bqm.add_variable(c2i[5, 3, 0, 3], 0)

bqm.add_variable(c2i[5, 3, 1, 0], 1)
bqm.add_variable(c2i[5, 3, 1, 1], 0)
bqm.add_variable(c2i[5, 3, 1, 2], 0)
# bqm.add_variable(c2i[5, 3, 1, 3], 0)


# bqm.add_variable(c2i[5, 4, 0, 0], 0)
# bqm.add_variable(c2i[5, 4, 0, 1], 0)
# bqm.add_variable(c2i[5, 4, 0, 2], 0)
# bqm.add_variable(c2i[5, 4, 0, 3], 0)

# bqm.add_variable(c2i[5, 4, 1, 0], 0)
# bqm.add_variable(c2i[5, 4, 1, 1], 0)
# bqm.add_variable(c2i[5, 4, 1, 2], 0)
# bqm.add_variable(c2i[5, 4, 1, 3], 0)




bqm.add_interaction(c2i[0, 0, 0, 0], c2i[0, 0, 1, 0], 1)
bqm.add_interaction(c2i[0, 0, 0, 0], c2i[0, 0, 1, 2], -1)
bqm.add_interaction(c2i[0, 0, 0, 0], c2i[1, 0, 0, 0], -1)
bqm.add_interaction(c2i[0, 0, 0, 1], c2i[0, 0, 1, 0], 1)
bqm.add_interaction(c2i[0, 0, 0, 1], c2i[0, 0, 1, 2], 0.5)
bqm.add_interaction(c2i[0, 0, 0, 1], c2i[1, 0, 0, 1], -1)
bqm.add_interaction(c2i[0, 0, 1, 2], c2i[0, 1, 1, 2], -1)
bqm.add_interaction(c2i[0, 1, 1, 2], c2i[0, 2, 1, 2], -1)
bqm.add_interaction(c2i[0, 2, 0, 2], c2i[0, 2, 1, 0], 1)
bqm.add_interaction(c2i[0, 2, 0, 2], c2i[0, 2, 1, 2], -1)
bqm.add_interaction(c2i[0, 2, 0, 2], c2i[1, 2, 0, 2], -1)
bqm.add_interaction(c2i[0, 2, 0, 3], c2i[0, 2, 1, 0], 1)
bqm.add_interaction(c2i[0, 2, 0, 3], c2i[0, 2, 1, 2], 0.5)
bqm.add_interaction(c2i[0, 2, 0, 3], c2i[1, 2, 0, 3], -1)
bqm.add_interaction(c2i[0, 2, 1, 2], c2i[0, 3, 1, 2], -1)
bqm.add_interaction(c2i[0, 3, 1, 2], c2i[0, 4, 1, 2], -1)
bqm.add_interaction(c2i[0, 4, 0, 0], c2i[0, 4, 1, 0], 1)
bqm.add_interaction(c2i[0, 4, 0, 0], c2i[0, 4, 1, 1], -1)
bqm.add_interaction(c2i[0, 4, 0, 0], c2i[0, 4, 1, 2], -1)
# bqm.add_interaction(c2i[0, 4, 0, 0], c2i[1, 4, 0, 0], -1)
bqm.add_interaction(c2i[0, 4, 0, 1], c2i[0, 4, 1, 0], 1)
# bqm.add_interaction(c2i[0, 4, 0, 1], c2i[0, 4, 1, 1], -1)
bqm.add_interaction(c2i[0, 4, 0, 1], c2i[0, 4, 1, 2], 0.5)
bqm.add_interaction(c2i[0, 4, 0, 1], c2i[1, 4, 0, 1], -1)
bqm.add_interaction(c2i[1, 0, 0, 0], c2i[1, 0, 1, 2], -1)
# bqm.add_interaction(c2i[1, 0, 0, 0], c2i[2, 0, 0, 0], -1)
# bqm.add_interaction(c2i[1, 0, 0, 1], c2i[1, 0, 1, 2], -1)
bqm.add_interaction(c2i[1, 0, 0, 1], c2i[2, 0, 0, 1], -1)
bqm.add_interaction(c2i[1, 0, 1, 2], c2i[1, 1, 1, 2], -1)
bqm.add_interaction(c2i[1, 1, 0, 0], c2i[1, 1, 1, 0], 1)
bqm.add_interaction(c2i[1, 1, 0, 0], c2i[1, 1, 1, 1], 1)
bqm.add_interaction(c2i[1, 1, 0, 0], c2i[1, 1, 1, 2], 1)
bqm.add_interaction(c2i[1, 1, 0, 1], c2i[1, 1, 1, 0], -1)
# bqm.add_interaction(c2i[1, 1, 0, 1], c2i[1, 1, 1, 1], -1)
bqm.add_interaction(c2i[1, 1, 0, 1], c2i[1, 1, 1, 2], 1)
bqm.add_interaction(c2i[1, 1, 0, 2], c2i[1, 1, 1, 0], -1)
bqm.add_interaction(c2i[1, 1, 0, 2], c2i[1, 1, 1, 1], 1)
bqm.add_interaction(c2i[1, 1, 0, 2], c2i[1, 1, 1, 2], -1)
bqm.add_interaction(c2i[1, 1, 0, 2], c2i[2, 1, 0, 2], -1)
bqm.add_interaction(c2i[1, 1, 1, 0], c2i[1, 2, 1, 0], -1)
bqm.add_interaction(c2i[1, 1, 1, 1], c2i[1, 2, 1, 1], -1)
# bqm.add_interaction(c2i[1, 1, 1, 2], c2i[1, 2, 1, 2], -1)
# bqm.add_interaction(c2i[1, 2, 0, 0], c2i[1, 2, 1, 0], -1)
bqm.add_interaction(c2i[1, 2, 0, 0], c2i[1, 2, 1, 1], -1)
# bqm.add_interaction(c2i[1, 2, 0, 0], c2i[1, 2, 1, 2], -1)
bqm.add_interaction(c2i[1, 2, 0, 0], c2i[2, 2, 0, 0], -1)
bqm.add_interaction(c2i[1, 2, 0, 1], c2i[1, 2, 1, 0], -1)
# bqm.add_interaction(c2i[1, 2, 0, 1], c2i[1, 2, 1, 1], -1)
# bqm.add_interaction(c2i[1, 2, 0, 1], c2i[1, 2, 1, 2], -1)
bqm.add_interaction(c2i[1, 2, 0, 1], c2i[2, 2, 0, 1], -1)
# bqm.add_interaction(c2i[1, 2, 0, 2], c2i[1, 2, 1, 0], -1)
# bqm.add_interaction(c2i[1, 2, 0, 2], c2i[1, 2, 1, 1], -1)
bqm.add_interaction(c2i[1, 2, 0, 2], c2i[1, 2, 1, 2], -1)
# bqm.add_interaction(c2i[1, 2, 0, 3], c2i[1, 2, 1, 0], -1)
# bqm.add_interaction(c2i[1, 2, 0, 3], c2i[1, 2, 1, 1], -1)
# bqm.add_interaction(c2i[1, 2, 0, 3], c2i[1, 2, 1, 2], -1)
bqm.add_interaction(c2i[1, 2, 0, 3], c2i[2, 2, 0, 3], -1)
# bqm.add_interaction(c2i[1, 2, 1, 0], c2i[1, 3, 1, 0], -1)
# bqm.add_interaction(c2i[1, 2, 1, 1], c2i[1, 3, 1, 1], -1)
bqm.add_interaction(c2i[1, 2, 1, 2], c2i[1, 3, 1, 2], -1)
bqm.add_interaction(c2i[1, 3, 0, 0], c2i[1, 3, 1, 0], 1)
bqm.add_interaction(c2i[1, 3, 0, 0], c2i[1, 3, 1, 1], 1)
bqm.add_interaction(c2i[1, 3, 0, 0], c2i[1, 3, 1, 2], 1)
bqm.add_interaction(c2i[1, 3, 0, 1], c2i[1, 3, 1, 0], -1)
# bqm.add_interaction(c2i[1, 3, 0, 1], c2i[1, 3, 1, 1], -1)
bqm.add_interaction(c2i[1, 3, 0, 1], c2i[1, 3, 1, 2], 1)
bqm.add_interaction(c2i[1, 3, 0, 2], c2i[1, 3, 1, 0], -1)
bqm.add_interaction(c2i[1, 3, 0, 2], c2i[1, 3, 1, 1], 1)
bqm.add_interaction(c2i[1, 3, 0, 2], c2i[1, 3, 1, 2], -1)
bqm.add_interaction(c2i[1, 3, 0, 2], c2i[2, 3, 0, 2], -1)
bqm.add_interaction(c2i[1, 3, 1, 0], c2i[1, 4, 1, 0], -1)
bqm.add_interaction(c2i[1, 4, 0, 0], c2i[1, 4, 1, 0], -1)
bqm.add_interaction(c2i[1, 4, 0, 0], c2i[2, 4, 0, 0], -1)
# bqm.add_interaction(c2i[1, 4, 0, 1], c2i[1, 4, 1, 0], -1)
bqm.add_interaction(c2i[1, 4, 0, 1], c2i[2, 4, 0, 1], -1)
bqm.add_interaction(c2i[2, 0, 0, 0], c2i[2, 0, 1, 0], 1)
bqm.add_interaction(c2i[2, 0, 0, 0], c2i[2, 0, 1, 1], -1)
bqm.add_interaction(c2i[2, 0, 0, 0], c2i[3, 0, 0, 0], -1)
bqm.add_interaction(c2i[2, 0, 0, 1], c2i[2, 0, 1, 0], 1)
bqm.add_interaction(c2i[2, 0, 0, 1], c2i[2, 0, 1, 1], 0.5)
bqm.add_interaction(c2i[2, 0, 0, 1], c2i[3, 0, 0, 1], -1)
bqm.add_interaction(c2i[2, 0, 1, 1], c2i[2, 1, 1, 1], -1)
# bqm.add_interaction(c2i[2, 1, 0, 2], c2i[2, 1, 1, 1], -1)
bqm.add_interaction(c2i[2, 1, 0, 2], c2i[3, 1, 0, 2], -1)
bqm.add_interaction(c2i[2, 1, 1, 1], c2i[2, 2, 1, 1], -1)
# bqm.add_interaction(c2i[2, 2, 0, 0], c2i[2, 2, 1, 0], -1)
# bqm.add_interaction(c2i[2, 2, 0, 0], c2i[2, 2, 1, 1], -1)
bqm.add_interaction(c2i[2, 2, 0, 0], c2i[3, 2, 0, 0], -1)
bqm.add_interaction(c2i[2, 2, 0, 1], c2i[2, 2, 1, 0], 1)
bqm.add_interaction(c2i[2, 2, 0, 1], c2i[2, 2, 1, 1], -1)
# bqm.add_interaction(c2i[2, 2, 0, 1], c2i[3, 2, 0, 1], -1)
bqm.add_interaction(c2i[2, 2, 0, 3], c2i[2, 2, 1, 0], 1)
bqm.add_interaction(c2i[2, 2, 0, 3], c2i[2, 2, 1, 1], 0.5)
bqm.add_interaction(c2i[2, 2, 0, 3], c2i[3, 2, 0, 3], -1)
bqm.add_interaction(c2i[2, 2, 1, 1], c2i[2, 3, 1, 1], -1)
# bqm.add_interaction(c2i[2, 3, 0, 2], c2i[2, 3, 1, 1], -1)
bqm.add_interaction(c2i[2, 3, 0, 2], c2i[3, 3, 0, 2], -1)
bqm.add_interaction(c2i[2, 3, 1, 1], c2i[2, 4, 1, 1], -1)
bqm.add_interaction(c2i[2, 4, 0, 0], c2i[2, 4, 1, 0], 1)
bqm.add_interaction(c2i[2, 4, 0, 0], c2i[2, 4, 1, 1], -1)
# bqm.add_interaction(c2i[2, 4, 0, 0], c2i[3, 4, 0, 0], -1)
bqm.add_interaction(c2i[2, 4, 0, 1], c2i[2, 4, 1, 0], 1)
bqm.add_interaction(c2i[2, 4, 0, 1], c2i[2, 4, 1, 1], 0.5)
bqm.add_interaction(c2i[2, 4, 0, 1], c2i[3, 4, 0, 1], -1)
bqm.add_interaction(c2i[3, 0, 0, 0], c2i[3, 0, 1, 3], -1)
# bqm.add_interaction(c2i[3, 0, 0, 0], c2i[4, 0, 0, 0], -1)
# bqm.add_interaction(c2i[3, 0, 0, 1], c2i[3, 0, 1, 3], -1)
bqm.add_interaction(c2i[3, 0, 0, 1], c2i[4, 0, 0, 1], -1)
bqm.add_interaction(c2i[3, 0, 1, 3], c2i[3, 1, 1, 3], -1)
bqm.add_interaction(c2i[3, 1, 0, 0], c2i[3, 1, 1, 0], -1)
bqm.add_interaction(c2i[3, 1, 0, 0], c2i[3, 1, 1, 1], 1)
bqm.add_interaction(c2i[3, 1, 0, 0], c2i[3, 1, 1, 2], -1)
bqm.add_interaction(c2i[3, 1, 0, 0], c2i[3, 1, 1, 3], 1)
bqm.add_interaction(c2i[3, 1, 0, 1], c2i[3, 1, 1, 0], 1)
bqm.add_interaction(c2i[3, 1, 0, 1], c2i[3, 1, 1, 1], 1)
bqm.add_interaction(c2i[3, 1, 0, 1], c2i[3, 1, 1, 2], -1)
bqm.add_interaction(c2i[3, 1, 0, 1], c2i[3, 1, 1, 3], -1)
# bqm.add_interaction(c2i[3, 1, 0, 2], c2i[3, 1, 1, 0], -1)
# bqm.add_interaction(c2i[3, 1, 0, 2], c2i[3, 1, 1, 1], -1)
bqm.add_interaction(c2i[3, 1, 0, 2], c2i[3, 1, 1, 2], 1)
# bqm.add_interaction(c2i[3, 1, 0, 2], c2i[3, 1, 1, 3], -1)
bqm.add_interaction(c2i[3, 1, 0, 3], c2i[3, 1, 1, 0], -1)
bqm.add_interaction(c2i[3, 1, 0, 3], c2i[3, 1, 1, 1], 1)
bqm.add_interaction(c2i[3, 1, 0, 3], c2i[3, 1, 1, 2], 1)
bqm.add_interaction(c2i[3, 1, 0, 3], c2i[3, 1, 1, 3], -1)
bqm.add_interaction(c2i[3, 1, 0, 3], c2i[4, 1, 0, 3], -1)
bqm.add_interaction(c2i[3, 1, 1, 0], c2i[3, 2, 1, 0], -1)
bqm.add_interaction(c2i[3, 1, 1, 1], c2i[3, 2, 1, 1], -1)
# bqm.add_interaction(c2i[3, 1, 1, 3], c2i[3, 2, 1, 3], -1)
# bqm.add_interaction(c2i[3, 2, 0, 0], c2i[3, 2, 1, 0], -1)
# bqm.add_interaction(c2i[3, 2, 0, 0], c2i[3, 2, 1, 1], -1)
bqm.add_interaction(c2i[3, 2, 0, 0], c2i[3, 2, 1, 3], -1)
# bqm.add_interaction(c2i[3, 2, 0, 0], c2i[4, 2, 0, 0], -1)
# bqm.add_interaction(c2i[3, 2, 0, 1], c2i[3, 2, 1, 0], -1)
bqm.add_interaction(c2i[3, 2, 0, 1], c2i[3, 2, 1, 1], -1)
# bqm.add_interaction(c2i[3, 2, 0, 1], c2i[3, 2, 1, 3], -1)
bqm.add_interaction(c2i[3, 2, 0, 1], c2i[4, 2, 0, 1], -1)
bqm.add_interaction(c2i[3, 2, 0, 2], c2i[3, 2, 1, 0], -1)
# bqm.add_interaction(c2i[3, 2, 0, 2], c2i[3, 2, 1, 1], -1)
# bqm.add_interaction(c2i[3, 2, 0, 2], c2i[3, 2, 1, 3], -1)
bqm.add_interaction(c2i[3, 2, 0, 2], c2i[4, 2, 0, 2], -1)
# bqm.add_interaction(c2i[3, 2, 0, 3], c2i[3, 2, 1, 0], -1)
# bqm.add_interaction(c2i[3, 2, 0, 3], c2i[3, 2, 1, 1], -1)
# bqm.add_interaction(c2i[3, 2, 0, 3], c2i[3, 2, 1, 3], -1)
bqm.add_interaction(c2i[3, 2, 0, 3], c2i[4, 2, 0, 3], -1)
# bqm.add_interaction(c2i[3, 2, 1, 0], c2i[3, 3, 1, 0], -1)
# bqm.add_interaction(c2i[3, 2, 1, 1], c2i[3, 3, 1, 1], -1)
bqm.add_interaction(c2i[3, 2, 1, 3], c2i[3, 3, 1, 3], -1)
bqm.add_interaction(c2i[3, 3, 0, 0], c2i[3, 3, 1, 0], -1)
bqm.add_interaction(c2i[3, 3, 0, 0], c2i[3, 3, 1, 1], -1)
bqm.add_interaction(c2i[3, 3, 0, 0], c2i[3, 3, 1, 2], 1)
bqm.add_interaction(c2i[3, 3, 0, 0], c2i[3, 3, 1, 3], 1)
bqm.add_interaction(c2i[3, 3, 0, 1], c2i[3, 3, 1, 0], 1)
bqm.add_interaction(c2i[3, 3, 0, 1], c2i[3, 3, 1, 1], -1)
bqm.add_interaction(c2i[3, 3, 0, 1], c2i[3, 3, 1, 2], 1)
bqm.add_interaction(c2i[3, 3, 0, 1], c2i[3, 3, 1, 3], -1)
# bqm.add_interaction(c2i[3, 3, 0, 2], c2i[3, 3, 1, 0], -1)
bqm.add_interaction(c2i[3, 3, 0, 2], c2i[3, 3, 1, 1], 1)
# bqm.add_interaction(c2i[3, 3, 0, 2], c2i[3, 3, 1, 2], -1)
# bqm.add_interaction(c2i[3, 3, 0, 2], c2i[3, 3, 1, 3], -1)
# bqm.add_interaction(c2i[3, 3, 0, 2], c2i[4, 3, 0, 2], -1)
bqm.add_interaction(c2i[3, 3, 0, 3], c2i[3, 3, 1, 0], -1)
bqm.add_interaction(c2i[3, 3, 0, 3], c2i[3, 3, 1, 1], 1)
bqm.add_interaction(c2i[3, 3, 0, 3], c2i[3, 3, 1, 2], 1)
bqm.add_interaction(c2i[3, 3, 0, 3], c2i[3, 3, 1, 3], -1)
bqm.add_interaction(c2i[3, 3, 0, 3], c2i[4, 3, 0, 3], -1)
bqm.add_interaction(c2i[3, 3, 1, 0], c2i[3, 4, 1, 0], -1)
bqm.add_interaction(c2i[3, 4, 0, 0], c2i[3, 4, 1, 0], -1)
bqm.add_interaction(c2i[3, 4, 0, 0], c2i[4, 4, 0, 0], -1)
# bqm.add_interaction(c2i[3, 4, 0, 1], c2i[3, 4, 1, 0], -1)
bqm.add_interaction(c2i[3, 4, 0, 1], c2i[4, 4, 0, 1], -1)
bqm.add_interaction(c2i[4, 0, 0, 0], c2i[4, 0, 1, 0], 1)
bqm.add_interaction(c2i[4, 0, 0, 0], c2i[4, 0, 1, 1], -1)
bqm.add_interaction(c2i[4, 0, 0, 0], c2i[5, 0, 0, 0], -1)
bqm.add_interaction(c2i[4, 0, 0, 1], c2i[4, 0, 1, 0], 1)
bqm.add_interaction(c2i[4, 0, 0, 1], c2i[4, 0, 1, 1], 0.5)
bqm.add_interaction(c2i[4, 0, 1, 1], c2i[4, 1, 1, 1], -1)
# bqm.add_interaction(c2i[4, 1, 0, 3], c2i[4, 1, 1, 1], -1)
bqm.add_interaction(c2i[4, 1, 0, 3], c2i[5, 1, 0, 3], -1)
bqm.add_interaction(c2i[4, 1, 1, 1], c2i[4, 2, 1, 1], -1)
bqm.add_interaction(c2i[4, 2, 0, 0], c2i[4, 2, 1, 0], -1)
# bqm.add_interaction(c2i[4, 2, 0, 0], c2i[4, 2, 1, 1], -1)
# bqm.add_interaction(c2i[4, 2, 0, 0], c2i[4, 2, 1, 2], -1)
bqm.add_interaction(c2i[4, 2, 0, 0], c2i[5, 2, 0, 0], -1)
# bqm.add_interaction(c2i[4, 2, 0, 1], c2i[4, 2, 1, 0], -1)
# bqm.add_interaction(c2i[4, 2, 0, 1], c2i[4, 2, 1, 1], -1)
# bqm.add_interaction(c2i[4, 2, 0, 1], c2i[4, 2, 1, 2], -1)
bqm.add_interaction(c2i[4, 2, 0, 1], c2i[5, 2, 0, 1], -1)
# bqm.add_interaction(c2i[4, 2, 0, 2], c2i[4, 2, 1, 0], -1)
bqm.add_interaction(c2i[4, 2, 0, 2], c2i[4, 2, 1, 1], -1)
bqm.add_interaction(c2i[4, 2, 0, 2], c2i[4, 2, 1, 2], 1)
# bqm.add_interaction(c2i[4, 2, 0, 3], c2i[4, 2, 1, 0], -1)
bqm.add_interaction(c2i[4, 2, 0, 3], c2i[4, 2, 1, 1], 0.5)
bqm.add_interaction(c2i[4, 2, 0, 3], c2i[4, 2, 1, 2], 1)
bqm.add_interaction(c2i[4, 2, 1, 0], c2i[4, 3, 1, 0], -1)
bqm.add_interaction(c2i[4, 2, 1, 1], c2i[4, 3, 1, 1], -1)
bqm.add_interaction(c2i[4, 3, 0, 2], c2i[4, 3, 1, 0], -1)
# bqm.add_interaction(c2i[4, 3, 0, 2], c2i[4, 3, 1, 1], -1)
bqm.add_interaction(c2i[4, 3, 0, 2], c2i[5, 3, 0, 2], -1)
# bqm.add_interaction(c2i[4, 3, 0, 3], c2i[4, 3, 1, 0], -1)
# bqm.add_interaction(c2i[4, 3, 0, 3], c2i[4, 3, 1, 1], -1)
bqm.add_interaction(c2i[4, 3, 0, 3], c2i[5, 3, 0, 3], -1)
# bqm.add_interaction(c2i[4, 3, 1, 0], c2i[4, 4, 1, 0], -1)
bqm.add_interaction(c2i[4, 3, 1, 1], c2i[4, 4, 1, 1], -1)
bqm.add_interaction(c2i[4, 4, 0, 0], c2i[4, 4, 1, 0], 1)
bqm.add_interaction(c2i[4, 4, 0, 0], c2i[4, 4, 1, 1], -1)
bqm.add_interaction(c2i[4, 4, 0, 1], c2i[4, 4, 1, 0], 1)
bqm.add_interaction(c2i[4, 4, 0, 1], c2i[4, 4, 1, 1], 0.5)
bqm.add_interaction(c2i[5, 0, 0, 0], c2i[5, 0, 1, 3], -1)
bqm.add_interaction(c2i[5, 0, 1, 3], c2i[5, 1, 1, 3], -1)
bqm.add_interaction(c2i[5, 1, 0, 0], c2i[5, 1, 1, 0], -1)
bqm.add_interaction(c2i[5, 1, 0, 0], c2i[5, 1, 1, 1], -1)
bqm.add_interaction(c2i[5, 1, 0, 0], c2i[5, 1, 1, 2], 1)
bqm.add_interaction(c2i[5, 1, 0, 0], c2i[5, 1, 1, 3], 1)
bqm.add_interaction(c2i[5, 1, 0, 1], c2i[5, 1, 1, 0], 1)
bqm.add_interaction(c2i[5, 1, 0, 1], c2i[5, 1, 1, 1], -1)
bqm.add_interaction(c2i[5, 1, 0, 1], c2i[5, 1, 1, 2], 1)
bqm.add_interaction(c2i[5, 1, 0, 1], c2i[5, 1, 1, 3], -1)
bqm.add_interaction(c2i[5, 1, 0, 2], c2i[5, 1, 1, 0], -1)
bqm.add_interaction(c2i[5, 1, 0, 2], c2i[5, 1, 1, 1], 1)
bqm.add_interaction(c2i[5, 1, 0, 2], c2i[5, 1, 1, 2], 1)
bqm.add_interaction(c2i[5, 1, 0, 2], c2i[5, 1, 1, 3], -1)
# bqm.add_interaction(c2i[5, 1, 0, 3], c2i[5, 1, 1, 0], -1)
bqm.add_interaction(c2i[5, 1, 0, 3], c2i[5, 1, 1, 1], 1)
# bqm.add_interaction(c2i[5, 1, 0, 3], c2i[5, 1, 1, 2], -1)
# bqm.add_interaction(c2i[5, 1, 0, 3], c2i[5, 1, 1, 3], -1)
bqm.add_interaction(c2i[5, 1, 1, 0], c2i[5, 2, 1, 0], -1)
# bqm.add_interaction(c2i[5, 1, 1, 2], c2i[5, 2, 1, 2], -1)
bqm.add_interaction(c2i[5, 2, 0, 0], c2i[5, 2, 1, 0], -1)
# bqm.add_interaction(c2i[5, 2, 0, 0], c2i[5, 2, 1, 2], -1)
# bqm.add_interaction(c2i[5, 2, 0, 1], c2i[5, 2, 1, 0], -1)
bqm.add_interaction(c2i[5, 2, 0, 1], c2i[5, 2, 1, 2], -1)
# bqm.add_interaction(c2i[5, 2, 1, 0], c2i[5, 3, 1, 0], -1)
bqm.add_interaction(c2i[5, 2, 1, 2], c2i[5, 3, 1, 2], -1)
bqm.add_interaction(c2i[5, 3, 0, 0], c2i[5, 3, 1, 0], 1)
bqm.add_interaction(c2i[5, 3, 0, 0], c2i[5, 3, 1, 1], 1)
bqm.add_interaction(c2i[5, 3, 0, 0], c2i[5, 3, 1, 2], 1)
# bqm.add_interaction(c2i[5, 3, 0, 1], c2i[5, 3, 1, 0], -1)
bqm.add_interaction(c2i[5, 3, 0, 1], c2i[5, 3, 1, 1], 1)
bqm.add_interaction(c2i[5, 3, 0, 1], c2i[5, 3, 1, 2], -1)
bqm.add_interaction(c2i[5, 3, 0, 2], c2i[5, 3, 1, 0], 1)
bqm.add_interaction(c2i[5, 3, 0, 2], c2i[5, 3, 1, 1], -1)
bqm.add_interaction(c2i[5, 3, 0, 2], c2i[5, 3, 1, 2], -1)
# bqm.add_interaction(c2i[5, 3, 0, 3], c2i[5, 3, 1, 0], -1)
bqm.add_interaction(c2i[5, 3, 0, 3], c2i[5, 3, 1, 1], -1)
# bqm.add_interaction(c2i[5, 3, 0, 3], c2i[5, 3, 1, 2], -1)




embedding = {'a0':  {c2i[0, 4, 0, 1], c2i[1, 4, 0, 1], c2i[2, 4, 0, 1], c2i[3, 4, 0, 1], c2i[4, 4, 0, 1]},
             'a1':  {c2i[0, 2, 0, 3], c2i[1, 2, 0, 3], c2i[2, 2, 0, 3], c2i[3, 2, 0, 3], c2i[4, 2, 0, 3]},
             'a2':  {c2i[0, 0, 0, 1], c2i[1, 0, 0, 1], c2i[2, 0, 0, 1], c2i[3, 0, 0, 1], c2i[4, 0, 0, 1]},
             'b0':  {c2i[0, 0, 1, 2], c2i[0, 1, 1, 2], c2i[0, 2, 1, 2], c2i[0, 3, 1, 2], c2i[0, 4, 1, 2]},
             'b1':  {c2i[2, 0, 1, 1], c2i[2, 1, 1, 1], c2i[2, 2, 1, 1], c2i[2, 3, 1, 1], c2i[2, 4, 1, 1]},
             'b2':  {c2i[4, 0, 1, 1], c2i[4, 1, 1, 1], c2i[4, 2, 1, 1], c2i[4, 3, 1, 1], c2i[4, 4, 1, 1]},
             'c10': {c2i[1, 3, 0, 2], c2i[2, 3, 0, 2], c2i[3, 3, 0, 2]},
             'c11': {c2i[1, 1, 0, 2], c2i[2, 1, 0, 2], c2i[3, 1, 0, 2]},
             'c21': {c2i[3, 3, 0, 3], c2i[4, 3, 0, 3], c2i[5, 3, 0, 3], c2i[5, 3, 1, 1]},
             'c22': {c2i[3, 1, 0, 3], c2i[4, 1, 0, 3], c2i[5, 1, 0, 3]},
             'c32': {c2i[5, 1, 1, 0], c2i[5, 2, 1, 0], c2i[5, 2, 0, 0], c2i[4, 2, 0, 0],
                     c2i[4, 2, 1, 0], c2i[4, 3, 1, 0], c2i[4, 3, 0, 2], c2i[5, 3, 0, 2]},
             'p0':  {c2i[0, 4, 0, 0], c2i[0, 4, 1, 1]},
             'p1':  {c2i[1, 3, 1, 1]},
             'p2':  {c2i[3, 3, 1, 2]},
             'p3':  {c2i[5, 3, 1, 0]},
             'p4':  {c2i[5, 1, 1, 2]},
             'p5':  {c2i[5, 1, 0, 2]},
             's00': {c2i[0, 2, 0, 2], c2i[1, 2, 0, 2], c2i[1, 2, 1, 2], c2i[1, 3, 1, 2]},
             's01': {c2i[0, 0, 0, 0], c2i[1, 0, 0, 0], c2i[1, 0, 1, 2], c2i[1, 1, 1, 2]},
             's11': {c2i[1, 1, 1, 1], c2i[1, 2, 1, 1], c2i[1, 2, 0, 0], c2i[2, 2, 0, 0],
                     c2i[3, 2, 0, 0], c2i[3, 2, 1, 3], c2i[3, 3, 1, 3]},
             's12': {c2i[2, 0, 0, 0], c2i[3, 0, 0, 0], c2i[3, 0, 1, 3], c2i[3, 1, 1, 3]},
             's22': {c2i[3, 1, 1, 1], c2i[3, 2, 1, 1], c2i[3, 2, 0, 1], c2i[4, 2, 0, 1],
                     c2i[5, 2, 0, 1], c2i[5, 2, 1, 2], c2i[5, 3, 1, 2]},
             's23': {c2i[4, 0, 0, 0], c2i[5, 0, 0, 0], c2i[5, 0, 1, 3], c2i[5, 1, 1, 3]},
             't01': {c2i[1, 3, 1, 0], c2i[1, 4, 1, 0], c2i[1, 4, 0, 0], c2i[2, 4, 0, 0]},
             't02': {c2i[3, 3, 1, 0], c2i[3, 4, 1, 0], c2i[3, 4, 0, 0], c2i[4, 4, 0, 0]},
             't11': {c2i[1, 1, 1, 0], c2i[1, 2, 1, 0], c2i[1, 2, 0, 1], c2i[2, 2, 0, 1]},
             't12': {c2i[3, 1, 1, 0], c2i[3, 2, 1, 0], c2i[3, 2, 0, 2], c2i[4, 2, 0, 2]},
             'aux0': {c2i[0, 0, 1, 0]},
             'aux1': {c2i[0, 2, 1, 0]},
             'aux2': {c2i[0, 4, 1, 0]},
             'aux3': {c2i[1, 1, 0, 0]},
             'aux4': {c2i[1, 1, 0, 1]},
             'aux5': {c2i[1, 3, 0, 0]},
             'aux6': {c2i[1, 3, 0, 1]},
             'aux7': {c2i[2, 0, 1, 0]},
             'aux8': {c2i[2, 2, 1, 0]},
             'aux9': {c2i[2, 4, 1, 0]},
             'aux10': {c2i[3, 1, 0, 0]},
             'aux11': {c2i[3, 1, 0, 1]},
             'aux12': {c2i[3, 1, 1, 2]},
             'aux13': {c2i[3, 3, 0, 0]},
             'aux14': {c2i[3, 3, 0, 1]},
             'aux15': {c2i[3, 3, 1, 1]},
             'aux16': {c2i[4, 0, 1, 0]},
             'aux17': {c2i[4, 2, 1, 2]},
             'aux18': {c2i[4, 4, 1, 0]},
             'aux19': {c2i[5, 1, 0, 0]},
             'aux20': {c2i[5, 1, 0, 1]},
             'aux21': {c2i[5, 1, 1, 1]},
             'aux22': {c2i[5, 3, 0, 0]},
             'aux23': {c2i[5, 3, 0, 1]}}




P = 21

# fix product qubits
fixed_variables = {}
fixed_variables.update(zip(('p5', 'p4', 'p3', 'p2', 'p1', 'p0'), "{:06b}".format(P)))
fixed_variables = {var: 1 if x == '1' else -1 for (var, x) in fixed_variables.items()}
for p, value in fixed_variables.items():
    for var in embedding[p]:
        bqm.fix_variable(var, value)
    embedding.pop(p)




# G = dnx.chimera_graph(16)
# G = dnx.chimera_graph(16, node_list=bqm.linear, edge_list=bqm.quadratic)
# # we need the mapping from each node in the target to its source node
# reverse_embedding = {}
# for v, chain in embedding.items():
#     for u in chain:
#         if u in reverse_embedding:
#             raise ValueError("target node {} assigned to more than one source node".format(u))
#         reverse_embedding[u] = v
# dnx.draw_chimera(G, linear_biases=bqm.linear, quadratic_biases=bqm.quadratic, labels=reverse_embedding)
# plt.show()




response = sampler.sample(bqm, num_reads=1000)
samples = embutil.unembed_samples(response, embedding, chain_break_method=embutil.discard)
