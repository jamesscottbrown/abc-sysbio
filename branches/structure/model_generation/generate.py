import sys, copy
import numpy as np

from network_defs import *
from writers import *

nspecies = 4
species = ['X', 'Y']

nX = node('X', [ reaction(['X'],['X','X'], "par(1)*ma1(X)") ])
nY = node('Y', [ reaction(['Y'],[], "par(2)*ma1(Y)") ])

eXY = edge( 'X', 'Y',
            [ reaction(['X','Y'],['Y','Y'], "switchp(par(5))*par(3)*ma2(X,Y)") ],
            [ reaction(['X','Y'],['X','X'], "switchn(par(6))*par(4)*ma2(X,Y)") ]
          )

ma1 = law( 'ma1', ['A'], 'A')
ma2 = law( 'ma2', ['A','B'], "A*B")

nodes = [nX, nY]
edges = [eXY]
laws = [ma1, ma2]


motif = network(species, nodes, edges, laws) 

mjp_code_writer(motif,'test_mjp.cu')
ode_code_writer(motif,'test_ode.cu')


