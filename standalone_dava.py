import networkx as nx
from standalone_dava_graphs import load_dgraph_from_lists

from standalone_dava_graphs import States

from standalone_dava_graphs import DiseaseGraph
from standalone_dava_graphs import plot_graph, bake_ids, to_directed
import time
import igraph
import copy
import numpy as np
import sys
import math



class DAVA_intervention():
    def __init__(self, budget, plotting=False):
        self.budget = budget
        self.plotting = plotting

    def intervene(self, dgraph, fast=True):
        bake_ids(dgraph.g)

        dgraph_copy = dgraph.copy()

        start = time.time()
        if fast==True:
            merged_dgraph, infected_master = self.merge(dgraph_copy)

            if is_tree(merged_dgraph.g):
                print("It is tree!")
                vacc_indexes = self.davaTREE(merged_dgraph, infected_master)
            else:
                #print("is not a tree.")
                vacc_indexes = self.davaNORMAL(merged_dgraph, infected_master)
        else:
            vacc_indexes = []
            for budget_used in range(self.budget):
                self.budget=1
                merged_dgraph, infected_master = self.merge(dgraph_copy)

                if is_tree(merged_dgraph.g):
                    print("It is tree!")
                    vacc = self.davaTREE(merged_dgraph, infected_master)
                else:
                    #print("is not a tree.")
                    vacc = self.davaNORMAL(merged_dgraph, infected_master)
                                
                assert len(vacc) <= 1
                if len(vacc) == 0: # the dominator tree had no neighbors (no susc. node reachable from any infected node)
                    print("Vaccinating further makes no difference.", budget_used, "vaccines used.")
                    break
                vacc_indexes.extend(vacc)
                dgraph_copy = dgraph.copy()
                dgraph_copy.g.delete_vertices(vacc_indexes) # remove nodes that are vaccinated up-to-now and restart the algorithm

        for r in vacc_indexes:
            dgraph.g.vs()[r]["vaccinated"] = True  # only used for plotting

        end = time.time()
        print(f"Running dava took {end-start} seconds.")
        return vacc_indexes

    # Combines all infected nodes into one infected node
    # New node has connections to all nodes next to an infected in graph
    def merge(self, dgraph: DiseaseGraph):
        g = dgraph.g
        assert len(dgraph.getInfectedNodes())>0

        susceptibles_next_to_infecteds = {} # maps susceptible nodes that are next to infected nodes to the weights on the edges to their infected neighbors
        for i in dgraph.getInfectedNodes():
            neighbours = dgraph.getSusceptibleNeighbors(i)
            for n in neighbours:
                if n.index in susceptibles_next_to_infecteds:
                    susceptibles_next_to_infecteds[n.index].append(g[i, n]) # the value being appended here is the weight between them
                else:
                    susceptibles_next_to_infecteds[n.index] = [g[i, n], ]

        to_delete_ids = [v.index for v in g.vs if v['state'] == States.I] # stored now such that infected master will not be deleted

        infected_master = g.add_vertex(state=States.I, x=0.5, y=0) # x and y are just used for plotting


        # If an S node had multiple I neighbors, the weight of the new edge is adjusted accordingly
        for (index, weights) in susceptibles_next_to_infecteds.items():
            prob = weights[0]
            for val in weights[1:]:
                prob = prob + (1-prob)*val # see "MERGE algorithm in the paper

            g.add_edge(infected_master.index, index, weight=prob)

        # Delete all old infected nodes
        g.delete_vertices(to_delete_ids)

        assert(len(dgraph.getInfectedNodes()) == 1)

        return dgraph, dgraph.getInfectedNodes()[0]


    # This is called if the merged graph is NOT a tree.
    # The function constructs the dominator tree from the merged graph and then calls davaTREE with it.
    def davaNORMAL(self, merged_dgraph, infected_master):
        dom_tree_dgraph = merged_dgraph.copy()
        dom_tree_dgraph.g.to_directed(mutual=True)

        # dom_tree is a list where each value is the id of the immediately dominating vertex of that index
        dom_tree = dom_tree_dgraph.g.dominator(
            vid=infected_master.index, mode='OUT')

        # delete old edges, we will now construct a tree by adding new (weighted) edges.
        dom_tree_dgraph.g.es.delete()

        # We now compute the edge weights, i.e. what the DAVA paper calls "maximum propagation path probability" (pdf p. 13)
        # to do this, we leverage igraph's shortest path implementation. But since that MINIMIZES the SUM of edge weights
        # (we want to MAXIMIZE the PRODUCT of edge weights), we preprocess the weights with -log(x) and later extract the
        # path lengths we need by computing exp(-x). This works since exp(-(-log(a) + -log(b))) = a*b and similarly for longer paths.
        edge_weights_for_shortest_paths = list(map(lambda x: -math.log(x),  merged_dgraph.g.es()["weight"]))
        shortest_paths = merged_dgraph.g.shortest_paths_dijkstra(source=infected_master.index, weights = edge_weights_for_shortest_paths)
        shortest_paths = shortest_paths[0] # list has only one entry since we only entered one source vertex

        for (index, dominator_index) in enumerate(dom_tree):
            # second check is for NaN (unreachable nodes in dom tree)
            if dominator_index != -1 and dominator_index == dominator_index:
                if dominator_index == infected_master.index:
                    weight = math.exp(- shortest_paths[index])
                else: # the DAVA paper notes that in the case of a dom. tree, this can be computed from other max prop. path probabilities: (pdf p.13)
                    weight = math.exp(shortest_paths[dominator_index] - shortest_paths[index]) # = math.exp(- shortest_paths[index])/ math.exp(- shortest_paths[dominator_index])

                #print("Adding edge", index, dominator_index, weight)
                dom_tree_dgraph.g.add_edge(dominator_index, index, weight=weight) 

        if self.plotting:
            plot_graph(dom_tree_dgraph, 'dava_dom_tree.png', layout=dom_tree_dgraph.g.layout_reingold_tilford_circular(mode="ALL",root=infected_master.index),)

        return self.davaTREE(dom_tree_dgraph, infected_master)

    def davaTREE(self, tree_dgraph, infected_master):
        benefits = []
        #print("Infected master has",len(tree_dgraph.g.neighbors(infected_master, mode="OUT")),"neighbors in the dom tree.")
        for n in tree_dgraph.g.neighbors(infected_master, mode="OUT"):
            weight = tree_dgraph.g[infected_master, n]
            partial = self.davaTREE_calPartial(
                tree_dgraph.g, n, parent_node=infected_master.index)
            benefits.append((n,weight*partial))

        benefits.sort(key=lambda x:x[1], reverse=True)
        # this list holds the dom-tree-indexes of the k (budget) nodes with highest benefit scores
        domtree_node_indexes = list(map(lambda x: x[0], benefits[:self.budget]))
          
        # since we removed some vertices along the way (during merging), we need to convert these indexes to the corresponding indexes in the input graph 
        originalgraph_indexes = []

        for indx in domtree_node_indexes:
            originalgraph_indexes.append(tree_dgraph.g.vs()[indx]["baked_index"])
            
        return originalgraph_indexes

    def davaTREE_calPartial(self, tree, start_node, parent_node=None):
        neighbors = tree.neighbors(start_node, mode="OUT")
        partial = 1
        for n in neighbors:
            if n == parent_node:
                continue
            partial += tree[start_node, n] * self.davaTREE_calPartial(tree, start_node=n,parent_node=start_node)

        return partial


# For some reason absolutely beyond me, the igraph library doesn't provide such a function.... -.-'
def is_tree(ig):
    if ig.get_edgelist() == []:
        return True
    nx_graph = nx.Graph(ig.get_edgelist())
    return nx.is_tree(nx_graph)



# For plotting pycairo needs to be installed, see https://pycairo.readthedocs.io/en/latest/getting_started.html
def dava_intervention(edge_list, infected_list=[0], recovered_list=[], k=1, fast=True, plotting=False):
    dgraph = load_dgraph_from_lists(edge_list, infected_list, recovered_list)
    if plotting:
        plot_graph(dgraph, name="dava_before_vaccination.png")

    # perform intervention
    di = DAVA_intervention(budget=k, plotting=plotting)
    node_ids = di.intervene(dgraph, fast=fast)

    if plotting:
        plot_graph(dgraph, name="dava_after_vaccination.png")

    return node_ids
    

if __name__ == '__main__':
    print("Starting DAVA standalone..")

    if len(sys.argv) < 4:
        print("Usage: %s <edgelist-file> <infected-file> <recovered-file> k. Give 'empty' to skip files." % (sys.argv[0]))
        exit()


    # load content of files
    graph_filepath = sys.argv[1]
    if graph_filepath == 'empty':
        print("Error: You have to give an edge list!")    
        exit()

    with open(graph_filepath, 'r') as graph_file:
        edge_list = list(map(lambda x: x.split(), graph_file.read().splitlines()))
        if len(edge_list[0]) == 2: # no weights
            edge_list = list(map(lambda e: (int(e[0]), int(e[1])), edge_list))
        elif len(edge_list[0]) == 3: # weights
            edge_list = list(map(lambda e: (int(e[0]), int(e[1]), float(e[2])), edge_list))

    infected_filepath = sys.argv[2]
    if infected_filepath == 'empty':
        infecteds=[]
    else:
        with open(infected_filepath, 'r') as infecteds_file:
            infecteds = list(map(int, infecteds_file.read().split()))

    recovered_filepath = sys.argv[3]
    if recovered_filepath == 'empty':
        recovereds=[]
    else:
        with open(recovered_filepath, 'r') as recovereds_file:
            recovereds = list(map(int, recovereds_file.read().split()))

    budget = int(sys.argv[4])

    vaccinated = dava_intervention(edge_list, infecteds, recovereds, k=budget, plotting=False, fast=False)
    print("Vaccinated:", vaccinated)
    
    F = open('dava_output.txt', 'w')
    F.writelines([str(v) + "\n" for v in vaccinated])
    F.close()
