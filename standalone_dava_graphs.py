from igraph import Graph, drawing
import igraph

import numpy as np
import copy

from enum import Enum


class States(Enum):
    S = "S"
    I = "I"
    R = "R"

class DiseaseGraph():
    def __init__(self, ig):
        self.g = ig

    def copy(self):
        return DiseaseGraph(ig=self.g.copy())

    def getInfectedNodes(self):
        return self.g.vs.select(state_eq=States.I)

    def getSusceptibleNodes(self):
        return self.g.vs.select(state_eq=States.S)

    def getNumNodesInState(self, state):
        return len(self.g.vs.select(state_eq=state))

    def getTotalNumNodes(self):
        return len(self.g.vs)

    def finishedSimulation(self):
        return self.getNumNodesInState(States.I) == 0

    def infect(self, target_node):
        if isinstance(target_node,int):
            target_node = self.g.vs[target_node]
        assert target_node['state'] == States.S
        target_node['state'] = States.I

    def recover(self, target_node):
        if isinstance(target_node,int):
            target_node = self.g.vs[target_node]

        target_node['state'] = States.R

    def susceptible(self, target_node):
        if isinstance(target_node,int):
            target_node = self.g.vs[target_node]

        target_node['state'] = States.S

    def setAllSusceptible(self):
        self.g.vs["state"] = States.S
        self.g.vs["vaccinated"] = False

    def infectRandomSusceptibleNodes(self, dgraph, number=1):
        nodes = np.random.choice(self.getSusceptibleNodes(), size=number)
        for node in nodes:
            node['state'] = States.I

    def getSusceptibleNeighbors(self, node):
        neighbors = node.neighbors()
        neighbors_susceptible = []
        for n in neighbors:
            if n['state'] == States.S:
                neighbors_susceptible.append(n)

        return neighbors_susceptible

        
def bake_ids(ig):
    for v in ig.vs:
        v["baked_index"] = v.index

    return

def to_directed(g, prob_SI):
    g.to_directed(mutual=True)
    g.es["weight"] = prob_SI

    return g

def generate_GRG_DiseaseGraph(n, p):
    igraph = Graph.GRG(n, p)
    if not igraph.is_connected():
        igraph = igraph.clusters().giant()
        print("Warning: generated graph not connected. Using largest subcomponent of size:", len(
            igraph.vs()))

    dgraph = DiseaseGraph(igraph)
    dgraph.setAllSusceptible()

    return dgraph


#Plotting
def prepare_graph_for_plotting(dgraph):
    dgraph.g.vs['color'] = 'black'
    for i in dgraph.g.vs():
        if i['vaccinated'] == True:
            i['color'] = 'red'
            continue

        if i['state'] == States.S:
            i['color'] = 'blue'
        if i['state'] == States.I:
            i['color'] = 'orange'
        if i['state'] == States.R:
            i['color'] = 'green'


# pycairo needs to be installed, see https://pycairo.readthedocs.io/en/latest/getting_started.html
def plot_graph(dgraph, name, bbox=(1000, 1000), **kwargs):
    prepare_graph_for_plotting(dgraph)
    drawing.plot(dgraph.g, target=name, bbox=bbox, **kwargs)


def load_dgraph_from_lists(edge_list, infected_list, recovered_list):
    DEFAULT_WEIGHT=0.2
    ig = igraph.Graph()

    if edge_list == []:
        return DiseaseGraph(ig)

    no_vertices = max( edge_list, key=lambda x: max(x[:2]))[1] # find highest number used in edge list
    ig.add_vertices(no_vertices+1)
    for e in edge_list:
        if len(e) == 2:    
            ig.add_edge(e[0],e[1],weight=DEFAULT_WEIGHT)
        elif len(e) == 3: 
            ig.add_edge(e[0],e[1],weight=e[2])
    
    dgraph = DiseaseGraph(ig)
    dgraph.setAllSusceptible()

    for inf in infected_list:
        dgraph.infect(inf)

    
    for rec in recovered_list:
        dgraph.recover(rec)

        # we cant just remove it since that would change the indexes of all nodes with higher indexes
        # thus we just delete all edges
        neighbors = dgraph.g.vs()[rec].neighbors()
        edges = []
        for n in neighbors:
            edges.append((rec,n.index))
            
        eids = dgraph.g.get_eids(pairs=edges)
        dgraph.g.delete_edges(eids)


    return dgraph