# Classes and methods for reactions, edges and nodes

import sys, copy
import numpy as np

class reaction:
    def __init__(self, reactants, products, rate):
        self.reactants = reactants
        self.products = products
        self.rate = rate
    
    def display(self):
        print self.reactants, '->', self.products
    def fdisplay(self,fout):
        print >>fout, self.reactants, '->', self.products

class edge:
    def __init__(self, start, end, reactionsp, reactionsn ):
        self.start = start
        self.end = end
        self.reactionsp = copy.deepcopy(reactionsp)
        self.reactionsn = copy.deepcopy(reactionsn)

    def display(self):
        print self.start,'->',self.end,':'
        for r in self.reactionsp:
            r.display()
        print self.start,'-|',self.end,':'
        for r in self.reactionsn:
            r.display()


class node:
    def __init__(self, name, reactions):
        self.name = name
        self.reactions = copy.deepcopy(reactions)
        
    def display(self):
        print self.name,':',
        for r in self.reactions:
            r.display()

class law:
    def __init__(self, name, inputs, output):
        self.name = name
        self.inputs = copy.deepcopy(inputs)
        self.output = copy.deepcopy(output)

## this helper function counts the instances of each species in the given list
def counter(species,list):
    ns = len(species)
    ret = np.zeros([ns])
    for i in range(ns):
        ret[i] = list.count(species[i])
    return(ret)

class network:
    def __init__(self, species, nodes, edges, laws, params):
        self.nodes = copy.deepcopy(nodes)
        self.edges = copy.deepcopy(edges)
        self.species = copy.deepcopy(species)
        self.laws = copy.deepcopy(laws)
        self.nspecies = len(self.species)
        self.params = params
        self.npar = len(params)

        # split the reactions into node reactions and edge reactions
        # since node reactions are always present
        self.total_reactions = []
        self.node_reactions = []
        self.edge_reactions = []
        self.edge_direction = []

        for i in self.nodes:
            for j in range(len(i.reactions)):
                self.node_reactions.append( i.reactions[j] )
                self.total_reactions.append( i.reactions[j] )
        for i in edges:
            for j in range(len(i.reactionsp)):
                self.edge_reactions.append( i.reactionsp[j] )
                self.edge_direction.append(1)
                self.total_reactions.append( i.reactionsp[j] )
            for j in range(len(i.reactionsn)):
                self.edge_reactions.append( i.reactionsn[j] )
                self.edge_direction.append(-1)
                self.total_reactions.append( i.reactionsn[j] )

        self.nnode_reactions = len(self.node_reactions)
        self.nedge_reactions = len(self.edge_reactions)
        self.nreactions = len(self.total_reactions)

        ## calculate the stoichiometry matrix
        self.smatrix = np.zeros( [self.nreactions,self.nspecies] )
        self.rmatrix = np.zeros( [self.nreactions,self.nspecies], dtype=np.int32 )

        for i in range(self.nreactions):
            # print 'reaction'
            # total_reactions[i].display()
            u = counter(self.species,self.total_reactions[i].reactants)
            v = counter(self.species,self.total_reactions[i].products)
            # print 'species  :', species
            # print 'reactants:', u
            # print 'products :', v
            # print 'net      :', v-u
            # print np.shape(smatrix[i,:]), np.shape(v-u)
            self.rmatrix[i,:] = u
            self.smatrix[i,:] = v-u

        self.sdim = np.shape(self.smatrix)     
