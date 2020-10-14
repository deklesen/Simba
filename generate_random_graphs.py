import os
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob
import random


########################################################
# Random graph models
########################################################


# nice graph for visuals
def geom_graph(node_num=350, radius = 0.07):
    seed = 42
    for _ in range(10000):
        G = nx.random_geometric_graph(node_num, radius, seed=seed)
        node_num = G.number_of_nodes()
        pos = nx.get_node_attributes(G, 'pos')
        node_pos = [(pos[i][0], pos[i][1]) for i in range(node_num)]
        G = nx.convert_node_labels_to_integers(G)
        edges = [e for e in G.edges]
        #print("Size of largest connected subgraph:",len(max(nx.connected_components(G), key=len)))
        if nx.is_connected(G):
            return G
        radius *= 1.1
        seed += 1
    print('failed graph generation geom_graph')
    return None


def fuzzy_geom_graph(size, radius, deg, ret_coords=True, force_connected=True):
    for _ in range(10000):
        # sample coordinates
        x, y = coords = np.random.rand(2, size) / radius

        # build the adjacency matrix
        adj = np.zeros((size, size)).astype(np.bool)
        for i, (xi, yi, di) in enumerate(zip(x, y, deg)):
            # sample neighbors based on euclidian distance
            p = np.exp(-np.sqrt((xi - x) ** 2 + (yi - y) ** 2))
            other_nodes = [k for k in range(size) if k != i]
            p = p[other_nodes]
            p /= p.sum()
            neighbors = np.random.choice(other_nodes, size=di, replace=False, p=p)
            adj[i, neighbors] = True
        adj |= adj.T

        G = nx.from_numpy_array(adj)
        if not force_connected or nx.is_connected(G):
            return G
    print('failed graph generation fuzzy_geom_graph')
    return None


def power_law_degrees(num_nodes=100, gamma=3.0, min_truncate=2, max_truncate=None):
    degree_distribution = [0] + [(k + 1) ** (-gamma) for k in range(num_nodes)]
    degree_distribution[:min_truncate] = [0.0] * min_truncate
    if max_truncate is not None:
        # max truncate and everything larger is zero
        degree_distribution[max_truncate:] = [0.0] * (len(degree_distribution) - max_truncate)
    assert (len(degree_distribution) == num_nodes + 1)
    z = np.sum(degree_distribution)
    degree_distribution = [p / z for p in degree_distribution]
    while True:
        degee_sequence = [np.random.choice(range(num_nodes + 1), p=degree_distribution) for _ in range(num_nodes)]
        if np.sum(degee_sequence) % 2 == 0:
            break
    return degee_sequence

def geom_powerlaw(num_nodes=100, gamma=2.2, min_truncate=2, max_truncate=None, radius=1.0):
    deg = power_law_degrees(num_nodes=num_nodes, gamma=gamma, min_truncate=min_truncate, max_truncate=max_truncate)
    G, coords = fuzzy_geom_graph(num_nodes, radius=radius, deg=deg, ret_coords=True, force_connected=True)
    return G, coords

def householdsuper_graph(num_nodes=100, hh_size_sample = None, gamma=1.8, min_truncate=2, max_truncate=15, seed=100):
    #assert(num_nodes % hh_size == 0)
    #hh_num = int(num_nodes / hh_size)
    #pl_graph, _ = power_law_graph(num_nodes=num_nodes, gamma=gamma, min_truncate=min_truncate, max_truncate=max_truncate, seed=seed)

    pl_graph = barabasi(num_nodes, m=3)

    nodes = list(pl_graph.nodes())

    def hh_size_sample_baseline():
        return random.choice([1,2,3,4,5])
    if hh_size_sample is None:
        hh_size_sample = hh_size_sample_baseline

    households = [list()]
    for n in nodes:
        h = households[-1]
        hh_size = hh_size_sample()
        if len(h) < hh_size:
            households[-1].append(n)
        else:
            households.append([n])

    for h in households:
        for n1 in h:
            for n2 in h:
                if n1>n2:
                    pl_graph.add_edge(n1,n2)

    hh_graph = nx.convert_node_labels_to_integers(pl_graph)

    return hh_graph




def household_graph(num_nodes=100, hh_size = 4, gamma=2.0, min_truncate=3, max_truncate=None, seed=100):
    # Important: seed only for power law
    assert(num_nodes % hh_size == 0)
    hh_num = int(num_nodes/hh_size)
    pl_graph, _ = power_law_graph(num_nodes=hh_num, gamma=gamma, min_truncate=min_truncate, max_truncate=max_truncate, seed=seed)

    pl_edges = list(pl_graph.edges())
    hh_edges = list()

    for node in pl_graph.nodes():
        new_nodes = [str(node)+'_'+str(ind) for ind in range(hh_size)]
        for n1 in new_nodes:
            for n2 in new_nodes:
                if n1 != n2:
                    hh_edges.append((n1,n2))

    for n1, n2 in pl_edges:
        n1_new = str(n1)+'_'+str(random.choice(range(hh_size)))
        n2_new = str(n2)+'_'+str(random.choice(range(hh_size)))
        hh_edges.append((n1_new, n2_new))

    hh_graph = nx.Graph(hh_edges)
    hh_graph = nx.convert_node_labels_to_integers(hh_graph)

    return hh_graph







# seed does not really work
def power_law_graph(num_nodes=100, gamma=2.2, min_truncate=2, max_truncate=None, seed=42):
    degree_distribution = [0] + [(k + 1) ** (-gamma) for k in range(num_nodes)]
    degree_distribution[:min_truncate] = [0.0] * min_truncate
    if max_truncate is not None:
        # max truncate and everything larger is zero
        degree_distribution[max_truncate:] = [0.0] * (len(degree_distribution) - max_truncate)
    assert (len(degree_distribution) == num_nodes + 1)
    z = np.sum(degree_distribution)
    degree_distribution = [p / z for p in degree_distribution]
    while True:
        seed += 1
        np.random.seed(seed)
        degee_sequence = [np.random.choice(range(num_nodes + 1), p=degree_distribution) for _ in range(num_nodes)]
        degee_sequence = [1 if d==0 else d for d in degee_sequence]
        if np.sum(degee_sequence) % 2 == 0:
            break

    for _ in range(100000):
        seed += 42
        contact_network = nx.configuration_model(degee_sequence, create_using=nx.Graph, seed=seed)
        contact_network.remove_edges_from(nx.selfloop_edges(contact_network))
        if nx.is_connected(contact_network):
            return contact_network

    print('Failed graph generation')
    return None


def geometric_configuration_model(degree_sequence, coordinates=None, ignore_location=False, mcmc=True):
    # if ignore_location then returns simple configuration model

    import itertools
    print('run geometric_configuration_model (might take a while) ', end='')
    for _ in range(10000):  # try with different coordinates
        print('.')
        G = nx.Graph()
        num_nodes = len(degree_sequence)
        if coordinates is None:
            coordinates = np.random.rand(2, num_nodes)

        stubs = [[v_i] * degree for v_i, degree in enumerate(degree_sequence)]
        stubs = list(itertools.chain.from_iterable(stubs))
        if len(stubs) % 2 != 0:
            raise ValueError("Sum of degree sequence must be even.")

        random.shuffle(stubs)
        for v_i in range(num_nodes):
            G.add_node(v_i)
            G.nodes[v_i]['pos'] = coordinates[:, v_i]

        for _ in range(len(stubs) * 2):
            if len(stubs) == 0:
                break
            v1, v2 = stubs[:2]
            if v1 == v2:
                random.shuffle(stubs)
            elif G.has_edge(v1, v2):
                random.shuffle(stubs)
            else:
                G.add_edge(v1, v2)
                stubs = stubs[2:]

        max_steps = len(G.edges) * 10000
        for i_step in range(max_steps):
            edges = list(G.edges())
            e1 = random.choice(edges)
            e2 = random.choice(edges)
            if len(set(list(e1) + list(e2))) == 4:  # make sure they do not share nodes
                e1_list = list(e1)
                e2_list = list(e2)
                random.shuffle(
                    e1_list)  # more suble trick to avoid bias (dont rewire primarily lower nodes and higher nodes with each other)
                random.shuffle(e2_list)
                new_edge1 = (e1_list[0], e2_list[0])
                new_edge2 = (e1_list[1], e2_list[1])

                if G.has_edge(*new_edge1) or G.has_edge(*new_edge2):
                    continue

                v1 = e1_list[0]
                v2 = e1_list[1]
                v3 = e2_list[0]
                v4 = e2_list[1]

                v1_pos = coordinates[:, v1]
                v2_pos = coordinates[:, v2]
                v3_pos = coordinates[:, v3]
                v4_pos = coordinates[:, v4]

                cost = np.linalg.norm(v1_pos - v2_pos) + np.linalg.norm(v3_pos - v4_pos)
                cost_new = np.linalg.norm(v1_pos - v3_pos) + np.linalg.norm(v2_pos - v4_pos)

                cost = cost ** 5
                cost_new = cost_new ** 5

                ratio = cost_new / (cost + cost_new)

                # if mcmc use probabilistic optimization
                if mcmc and ratio > random.random():
                    continue

                # if not mcmc simply compare old and new costs
                # if ignore_location, no rejection steps are performed
                if (not mcmc) and (not ignore_location) and cost_new > cost:
                    continue

                G.add_edge(*new_edge1)
                G.add_edge(*new_edge2)
                G.remove_edge(*e1)
                G.remove_edge(*e2)

                if (i_step > max_steps / 2):
                    print(nx.is_connected(G), stubs)

                if len(stubs) > 0:
                    v1, v2 = stubs[:2]
                    if v1 != v2 and not G.has_edge(v1, v2):
                        G.add_edge(v1, v2)
                        stubs = stubs[2:]

                if (i_step > max_steps / 100) and nx.is_connected(G) and len(
                        stubs) == 0:  # note that we need burn in period
                    return G
                if (i_step > max_steps / 20) and nx.is_connected(G):
                    return G
    print('Generation failed geometric_configuration_model')
    return None


def geometric_configuration_pl(num_nodes=100, gamma=2.2, min_truncate=2, max_truncate=None, mcmc=False):
    deg = power_law_degrees(num_nodes=num_nodes, gamma=gamma, min_truncate=min_truncate, max_truncate=max_truncate)
    G = geometric_configuration_model(degree_sequence=deg, mcmc=mcmc)
    return G

def geometric_configuration_regular(num_nodes=100, d=3, mcmc=False):
    deg = [d]*num_nodes
    G = geometric_configuration_model(degree_sequence=deg, mcmc=mcmc)
    return G



def test_household_graph():
    plt.clf()
    G, _ = householdsuper_graph(100,4)
    pos = nx.spring_layout(G)
    nx.draw_networkx_edges(G, pos=pos, width=0.5)
    nx.draw_networkx_nodes(G, pos=pos, alpha=0.5, node_size=40)
    plt.savefig('household_graph.pdf')

def test_geometric_configuration_model():
    plt.clf()
    G, _ = geometric_configuration_model([3] * 500, mcmc=False)
    pos = nx.spring_layout(G)
    pos = {v_i: G.nodes[v_i]['pos'] for v_i in G.nodes()}
    nx.draw_networkx_edges(G, pos=pos, width=0.5)
    nx.draw_networkx_nodes(G, pos=pos, alpha=0.5, node_size=40)
    plt.savefig('example_geom_config.pdf')


    plt.clf()
    G, _ = geometric_configuration_pl()
    pos = {v_i: G.nodes[v_i]['pos'] for v_i in G.nodes()}
    nx.draw_networkx_edges(G, pos=pos, width=0.5)
    nx.draw_networkx_nodes(G, pos=pos, alpha=0.5, node_size=40)
    plt.savefig('example_geom_config_pl.pdf')

########################################################
# From other file
########################################################


def geom_graph2(num_nodes=100, radius=0.1):
    seed = 42
    for _ in range(100000):
        G = nx.random_geometric_graph(num_nodes, radius, seed=seed)
        node_num = G.number_of_nodes()
        pos = nx.get_node_attributes(G, 'pos')
        node_pos = [(pos[i][0], pos[i][1]) for i in range(node_num)]
        if nx.is_connected(G):
            G = nx.convert_node_labels_to_integers(G)
            return G
        radius *= 1.1
        seed += 1
    print('failed graph generation geom_graph2')


def erdos_renyi(num_nodes=100, connect_prob=0.05, mean_degree=None):
    if mean_degree is not None:
        assert (connect_prob is None)
        connect_prob = mean_degree / num_nodes
    for seed in range(1000):
        seed += 333
        G = nx.fast_gnp_random_graph(n=num_nodes, p=connect_prob, seed=seed)
        if nx.is_connected(G):
            #print("Generated graph mean degree:",np.mean([x[1] for x in G.degree()]))
            G = nx.convert_node_labels_to_integers(G)
            return G
        connect_prob=connect_prob*1.05
    print('failed graph generation erdos_renyi')


def barabasi(num_nodes=100, m=2):
    for seed in range(10000):
        seed += 22
        G = nx.barabasi_albert_graph(n=num_nodes, m=m)
        if nx.is_connected(G):
            G = nx.convert_node_labels_to_integers(G)
            return G
        m=m*1.05
    print('failed graph generation')
    return None


def regular(num_nodes=100, d=5):
    G = nx.random_regular_graph(d=d, n=num_nodes)
    G = nx.convert_node_labels_to_integers(G)
    return G


def complete(num_nodes=100):
    G = nx.complete_graph(num_nodes)
    return G


def newman(num_nodes=100, k=6, p=0.2):
    for seed in range(10000):
        seed += 222
        G = nx.watts_strogatz_graph(n=num_nodes, k=k, p=p, seed=seed)
        if nx.is_connected(G):
            G = nx.convert_node_labels_to_integers(G)
            return G
    print('failed graph generation')
    return None


def grid_2d(dim=10):
    num_nodes = dim * dim
    G = nx.grid_2d_graph(n=dim, m=dim)
    G = nx.convert_node_labels_to_integers(G)
    return G
