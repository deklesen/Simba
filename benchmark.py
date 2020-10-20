from baselines import *
from generate_random_graphs import *
import timeit
from functools import partial
from pathlib import Path
import os, time
import pandas as pd

run_nr=1337#random.randint(0,100000)
print("Run number:",run_nr)


def call_rust(CG, init_infected, path, init_recovered = None, num_run=1000, infection_rate=2.0):
    exec = '/rust_code/rust_reject/target/release/rust_reject'
    path=path+'{type}.{fileformat}'
    graph = path.format(type='init_graph', fileformat='txt')
    out_traj = path.format(type='out_trajectory', fileformat='txt')
    out_tg = path.format(type='TG', fileformat='edgelist')


    nodes = sorted(list(CG.nodes()))
    assert(nodes[-1] == CG.number_of_nodes()-1)
    if init_recovered is not None:
        for n in init_recovered:
            assert(n not in init_infected)
        for n in init_infected:
            assert(n not in init_recovered)

    with open(graph, 'w') as f:
        for n in nodes:
            l = 'S'
            if n in init_infected: l = 'I'
            if init_recovered is not None and n in init_recovered: l = 'R'
            neigh_list = sorted(list(CG.neighbors(n)))
            neigh_list = ','.join([str(v) for v in neigh_list])
            f.write('{};{};{}\n'.format(n,l,neigh_list))

    #print('start rust')
    os.system('.{} {} {} {} {} {} > /dev/null'.format(exec, graph, out_traj, out_tg, infection_rate, num_run))
    #time.sleep(0.05)
    #print('rust ended, back in python')

    #traj_data = pd.read_csv(out_traj, sep=',')
    #TG = nx.read_weighted_edgelist(out_tg, create_using=nx.DiGraph(), delimiter=' ', nodetype=int)

    mean_values = np.loadtxt(out_traj+'.score')
    mean = np.mean(mean_values)
    stddev = np.std(mean_values)

    #number_of_times_infected = np.loadtxt(out_tg+'.intensity')
    #intensity = number_of_times_infected/num_run
    #for n in nodes:
    #    TG.add_node(n)
    #    TG.nodes[n]['intensity'] = intensity[n]
    #TG.nodes[len(intensity)]['intensity'] = 0 #dummy node

    # load eq
    #eq_dist = np.loadtxt(out_tg+'.solution')

    # plot graphs
    if False and CG.number_of_nodes() < PLOTLIM:
        #pass
        pos_cg = nx.kamada_kawai_layout(CG)
        pos_cg = [(pos_cg[i][0], pos_cg[i][1]) for i in range(CG.number_of_nodes())]
        pos_tg = pos_cg + [(0, 1.5 * np.max([p[1] for p in pos_cg]))]
        plot_contact_graph(CG, path.format(type='init_graph_'+type, fileformat='pdf'), init_infected, init_recovered=init_recovered, pos=pos_cg)
        plot_TG(TG, path.format(type='transm_graph_'+type, fileformat='pdf'), init_infected, eq_dist, number_of_times_infected, init_recovered=init_recovered, pos=pos_tg,  CG=CG, pos_cg=pos_cg)

    #print('output done in call rust')
    #return CG, TG, traj_data, final_unaffected_mean, eq_dist
    return mean, stddev




################################################
################################################
################################################
################################################
################################################
################################################
################################################


# Filter out None-graphs
# Transform graphs to largest subgraph
def filter_graphs(input):
    to_largest_subgraph = lambda graph: graph.subgraph(max(nx.connected_components(graph), key=len))
    a = dict(filter(lambda x: x[1] is not None, input.items()))
    return { name: to_largest_subgraph(g) for name, g in a.items()}

import numpy as np

def get_graph_data(graph):
    #print(graph)
    degrees = np.array(list(map(lambda x: x[1], graph.degree())))
    return {
        'num_nodes': graph.number_of_nodes(),
        'num_edges': graph.number_of_edges(),
        'density': nx.density(graph),
        'degree_mean': degrees.mean(),
        'degree_stddev': degrees.std(),
        'radius': nx.radius(graph),
        'diameter': nx.diameter(graph),
        'transitivity': nx.transitivity(graph),
        'average_clustering': nx.average_clustering(graph),
        'average_shortest_path_length': nx.average_shortest_path_length(graph), ##slow
    }

baselines={
    'random': baseline_random,
    #'dava': baseline_dava,
    'davaFast': baseline_davaf,
    'degree': baseline_degree,
    'betweenness': baseline_betweenness,
    'EVcentrality': baseline_EVcentrality,
    'closeness': baseline_closenessCentrality,
    'PageRank': baseline_Pagerank,
    'PersPageRank': baseline_PersPagerank
}

graphs=filter_graphs({
    **{f'geom_graph_{node_num}': geom_graph(node_num) for node_num in [150,300,500,1000,2000,2500]},
    **{f'householdsuper_{node_num}':householdsuper_graph(node_num) for node_num in [150,300,500,1000,2000,2500]},
    **{f'erdos_renyi_{node_num}': erdos_renyi(node_num) for node_num in [150,300,500,1000,2000,2500]},
    **{f'barabasi_{node_num}': barabasi(node_num) for node_num in [25,50,100,250,500,1000]},
    **{f'grid_2d{node_num}': grid_2d(node_num) for node_num in [10,20,30,50]},
    **{f'newman{node_num}': newman(node_num) for node_num in [10,20,30,50]},
    **{f'complete{node_num}': complete(node_num) for node_num in [10,20,30,50]},
    **{f'regular{node_num}': regular(node_num) for node_num in [10,20,30,50]},
})

infection_rates = [0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5]

init_infecteds_fraction=[0.025,0.05,0.075,0.1,0.15,0.25]
budgets_fraction=[0.025,0.05,0.075,0.1,0.15,0.25]
SIMULATION_RUNS = 2000

import pickle
import tqdm
import time
import sys

if __name__=='__main__':
    runner_id = int(sys.argv[1])-1
    total_runners = int(sys.argv[2])


    print(f"Found {len(baselines)} baselines and {len(graphs)} graphs.")
    num_experiments=len(graphs)*len(baselines)*len(budgets_fraction)*len(init_infecteds_fraction)*len(infection_rates)
    print("Doing", num_experiments, "experiments...")

    run_counter=1

    results = {}
    with tqdm.tqdm(total=num_experiments) as pbar:
        for infection_rate in infection_rates:
            for graph_name, graph in graphs.items():
                graph_data = {
                    **get_graph_data(graph),
                    'name': graph_name,
                }
                for iif in init_infecteds_fraction:
                    num_init_infected = int(len(graph)*iif)
                    init_infected = random.sample(list(range(len(graph))), num_init_infected)

                    for baseline_name, baseline_func in baselines.items():
                        for budget_fraction in budgets_fraction:   
                            run_counter+=1
                            pbar.update(1)
                            if not (run_counter % total_runners == runner_id):
                                continue


                            budget=int(len(graph)*budget_fraction)

                            path="output/benchmark/{run_nr}/infrate{infection_rate}/{graph_name}_inf{init_infected}/{baseline_name}_budget{budget}/"
                            outpath = path.format(run_nr=run_nr, infection_rate=infection_rate, graph_name=graph_name, init_infected=num_init_infected, baseline_name=baseline_name, budget=budget)
                            if not Path(outpath).exists():
                                Path(outpath).mkdir(parents=True)
                            
                            start_time = time.time()
                            vaccinated = baseline_func(graph.copy(), init_infected, budget)
                            duration = time.time() - start_time

                            score, stddev = call_rust(graph.copy(), init_infected, path=outpath, init_recovered=vaccinated, infection_rate=infection_rate, num_run=SIMULATION_RUNS)
                            
                            result_data = {'score_mean':score, 'score_stddev':stddev, 'duration':duration}

                            results[(infection_rate,graph_name,baseline_name,iif,budget)] = (graph_data,result_data)
                            


    with open("output/benchmark/{run_nr}/results_{runner_id}.pickle".format(run_nr=run_nr,runner_id=runner_id), 'wb') as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)