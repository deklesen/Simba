import networkx as nx
import logging
import numpy as np
import pandas as pd
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import time
import os

import generate_random_graphs as gg

PLOTLIM = 99
OPTIMSTEPS = 25


logging.basicConfig(level=logging.INFO, filename = time.strftime("simba_%Y-%m-%d.log"), format = '%(asctime)s  %(levelname)-10s %(message)s')
logging.getLogger().addHandler(logging.StreamHandler())
#logging.info("info")
(logging.getLogger()).handlers[0].flush()

def speed_test(num_node_list = None, degree_list = None):
    # Regular
    if num_node_list is None:
        num_node_list = [10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6]
    if degree_list is None:
        degree_list = [5, 20, 100]
    for num_nodes in num_node_list:
        for d in degree_list:
            if num_nodes == 100 and d == 100:
                num_nodes += 1
            init_infected=[0,1]
            CG=gg.regular(num_nodes=num_nodes, d=d)
            path = 'output/'+'regular_{}_{}_a'.format(num_nodes, d)+'/{type}.{fileformat}'
            os.system('mkdir output/')
            os.system('mkdir output/' + 'regular_{}_{}_a'.format(num_nodes, d))
            type = 'single_run'
            call_rust(CG, init_infected, path, type, init_recovered=None, num_run=1000, infection_rate=2.0)


            #analysis(gg.regular(num_nodes=num_nodes, d=d), 'regular_{}_{}_a'.format(num_nodes, d), infection_rate=3.0 / d, budget=2, init_infected=[0])
            #analysis(gg.regular(num_nodes=num_nodes, d=d), 'regular_{}_{}_b'.format(num_nodes, d), infection_rate=3.0 / d, budget=10, init_infected=[0, 1, 2, 3, 4])

def plot_optimization_summary(optimization_summary, outpath, no_vacc_score=None, random_score=None, davaf=None, dava=None, pagerank_baseline=None, pers_pagerank_baseline=None, degree_baseline=None):
    optimization_best = list()
    for i in range(len( optimization_summary['score'])):
        optimization_best.append(np.max( optimization_summary['score'][:i+1]))

    plt.clf()
    plt.plot(optimization_summary['step_i'],optimization_best, color=sns.xkcd_rgb['denim blue'], lw=2,
             alpha=0.5, zorder=80, label='Simba')


    max_value = np.max(optimization_summary['score'])
    
    #plt.plot([0, optimization_summary['step_i'][-1]], [max_value, max_value], color='white' ,lw=3,
    #         alpha=0.7, zorder=1)
    #plt.plot([0, optimization_summary['step_i'][-1]], [max_value, max_value], color=sns.xkcd_rgb['pinkish red'], lw=2,
    #         alpha=0.9, zorder=20, label='Simba (final)')

    plt.scatter(optimization_summary['step_i'], optimization_best, color=sns.xkcd_rgb['pinkish red'], edgecolors='white', zorder=100, linewidths=1.5, s=80, alpha=0.7)

    plt.ylim([0, 1.05])
    plt.xlim([0, optimization_summary['step_i'][-1]])
    plt.ylabel(r'$F(X)$', fontsize=26)
    plt.xlabel('Iteration', fontsize=26)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)

    plt.locator_params(axis='y', nbins=5)
    plt.locator_params(axis='x', nbins=5)


    #if no_vacc_score is not None:
    #    plt.plot([0, optimization_summary['step_i'][-1]], [no_vacc_score, no_vacc_score], color=sns.xkcd_rgb['green'], label='wo vaccination', alpha=0.7)


    if dava is not None:
        plt.plot([0, optimization_summary['step_i'][-1]], [dava, dava], color=sns.xkcd_rgb['orange'],
                 ls=':', label='DAVA', alpha=0.6, lw=6, zorder=4)

    if davaf is not None:
        plt.plot([0, optimization_summary['step_i'][-1]], [davaf, davaf], color=sns.xkcd_rgb['orange'],
                 ls=':', label='DAVA fast', alpha=0.3, lw=9, zorder=3)


    if pagerank_baseline is not None:
        plt.plot([0, optimization_summary['step_i'][-1]], [pagerank_baseline, pagerank_baseline], color=sns.xkcd_rgb['pinkish purple'],
                 ls='-.', label='PageRank', alpha=0.6, lw=5, zorder=4)

    if pers_pagerank_baseline is not None:
        plt.plot([0, optimization_summary['step_i'][-1]], [pers_pagerank_baseline, pers_pagerank_baseline], color=sns.xkcd_rgb['pinkish purple'],
                 ls='-.', label='Pers. PageRank', alpha=0.3, lw=7, zorder=3)


    if degree_baseline is not None:
        plt.plot([0, optimization_summary['step_i'][-1]], [degree_baseline, degree_baseline], color=sns.xkcd_rgb['dark blue'],
                 ls='--', label='Degree', alpha=0.8, lw=3, zorder=4)

    if random_score is not None:
        plt.plot([0, optimization_summary['step_i'][-1]], [random_score, random_score], color=sns.xkcd_rgb['aqua'],
                 ls='--', label='random', alpha=0.4, lw=7, zorder=3)

    plt.legend(fontsize='large', bbox_to_anchor=(1.04, 1), loc="upper left", frameon=False)
    plt.savefig(outpath, bbox_inches="tight")

    min1 = np.min(optimization_summary['score'])
    minv = np.min([min1, dava, davaf, pagerank_baseline, pers_pagerank_baseline, random_score, degree_baseline])
    max1 = np.max(optimization_summary['score'])
    maxv = np.max([max1, dava, davaf, pagerank_baseline, pers_pagerank_baseline, random_score, degree_baseline])

    plt.ylim([minv*0.975, maxv*1.02])
    plt.savefig(outpath.replace('.pdf', '_scaled.pdf'), bbox_inches="tight")

    #no_vacc_score=None, random_score=None, davaf=None, dava=None, pagerank_baseline=None, pers_pagerank_baseline=None, degree_baseline=None
    outpath_data = outpath.replace('.pdf', '_data.csv')
    data = {'no_vacc_score': [no_vacc_score], 'random_score':[random_score], 'davaf':[davaf], 'pagerank_baseline':[pagerank_baseline], 'pers_pagerank_baseline':[pers_pagerank_baseline], 'degree_baseline':[degree_baseline]}
    data_df = pd.DataFrame(data)
    data_df.to_csv(outpath_data)


def plot_contact_graph(CG, outpath, init_infected, init_recovered=None, pos=None):
    if init_recovered is None:
        init_recovered = list()
    if CG.number_of_nodes() > PLOTLIM:
        return
    plt.clf()
    if pos is None:
        pos = nx.kamada_kawai_layout(CG)

    node_pos_x = [pos[i][0] for  i in range(CG.number_of_nodes())]
    node_pos_y = [pos[i][1] for  i in range(CG.number_of_nodes())]
    node_pos_x_S = [pos[i][0] for  i in range(CG.number_of_nodes()) if i not in init_recovered and i not in init_infected]
    node_pos_y_S = [pos[i][1] for  i in range(CG.number_of_nodes()) if i not in init_recovered and i not in init_infected]
    node_pos_x_I = [pos[i][0] for  i in range(CG.number_of_nodes()) if i in init_infected]
    node_pos_y_I = [pos[i][1] for  i in range(CG.number_of_nodes()) if i in init_infected]
    node_pos_x_R = [pos[i][0] for  i in range(CG.number_of_nodes()) if i in init_recovered]
    node_pos_y_R = [pos[i][1] for  i in range(CG.number_of_nodes()) if i in init_recovered]

    plt.scatter(node_pos_x_S, node_pos_y_S, s=100, alpha=0.8, c='black',
                edgecolors='none', zorder=15)
    plt.scatter(node_pos_x_I, node_pos_y_I, s=100, alpha=0.8, c='red',
                edgecolors='none', zorder=15)
    plt.scatter(node_pos_x_R, node_pos_y_R, s=100, alpha=0.8, c='orange',
                edgecolors='none', zorder=15)

    # plot labels
    for i in range(len(node_pos_x)):
        plt.text(node_pos_x[i], node_pos_y[i], i, alpha=0.8, color='white', fontsize=4, horizontalalignment='center',
                 verticalalignment='center', zorder=20)


    lw = min(3, 1.0 / (len(CG.edges())+1) * 600)
    for e in CG.edges:
        pos_v1_x = node_pos_x[e[0]]
        pos_v1_y = node_pos_y[e[0]]
        pos_v2_x = node_pos_x[e[1]]
        pos_v2_y = node_pos_y[e[1]]
        plt.plot([pos_v1_x, pos_v2_x], [pos_v1_y, pos_v2_y], c='black', alpha=0.5, zorder=10, lw=lw)

    for p in ['top', 'right', 'bottom', 'left']:
        (plt.gca()).spines[p].set_visible(False)
    (plt.gca()).set_yticklabels([])
    (plt.gca()).set_xticklabels([])
    plt.xticks([], [])
    plt.yticks([], [])

    plt.savefig(outpath, bbox_inches='tight', dpi=300)
    node_pos = [(node_pos_x[i], node_pos_y[i]) for i in range(len(node_pos_x))]
    return node_pos


def plot_TG(TG, outpath, init_infected, eq_dist, number_of_times_infected, init_recovered=None, pos=None, CG=None, pos_cg=None):
    if init_recovered is None:
        init_recovered = list()
    if TG.number_of_nodes() > PLOTLIM:
        return
    plt.clf()
    if pos is None:
        pos = nx.kamada_kawai_layout(TG)

    eq_dist_filterd =  [eq_dist[i] for i in range(len(eq_dist)) if i not in init_infected+init_infected]
    options = {
        'node_color': eq_dist_filterd,
        'edge_color': ['black']*len(TG.edges()),
        'node_size': 70,
        'width': 1,
        'arrowstyle': '-|>',
        'arrowsize': 7,
        'alpha': 0.8,
        'with_labels': True,
        'edgecolors': 'none',
        "edge_cmap": plt.cm.YlGnBu,
    }
    nx.draw_networkx(TG, pos=pos, arrows=True, **options, nodelist=list())

    options['alpha'] = 0.8
    nx.draw_networkx(TG, pos=pos, arrows=True, **options, edgelist=list(), cmap=plt.cm.hot, nodelist = [i for i in range(len(eq_dist)) if i not in init_infected+init_infected])
    nx.draw_networkx(TG, pos=pos, arrows=True, edgelist=list(), node_color='orange',  nodelist = init_recovered, node_shape='^')
    nx.draw_networkx(TG, pos=pos, arrows=True, edgelist=list(), node_color='red', nodelist = init_infected, node_shape='^')

    # draw CG
    if CG is not None:
        nx.draw_networkx(CG, pos=pos_cg, nodelist=list(), alpha=0.5)

    for p in ['top', 'right', 'bottom', 'left']:
        (plt.gca()).spines[p].set_visible(False)
    plt.savefig(outpath, bbox_inches='tight', dpi=300)


import os, time
import pandas as pd
def call_rust(CG, init_infected, path, type, init_recovered = None, num_run=1000, infection_rate=2.0):
    exec = '/rust_code/rust_reject/target/release/rust_reject'
    graph = path.format(type='init_graph_'+type, fileformat='txt')
    out_traj = path.format(type='out_trajectory_'+type, fileformat='txt')
    out_tg = path.format(type='TG_'+type, fileformat='edgelist')

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

    print('start rust')
    os.system('.{} {} {} {} {} {}'.format(exec, graph, out_traj, out_tg, infection_rate, num_run))
    time.sleep(0.05)
    print('rust ended, back in python')

    traj_data = pd.read_csv(out_traj, sep=',')
    TG = nx.read_weighted_edgelist(out_tg, create_using=nx.DiGraph(), delimiter=' ', nodetype=int)

    mean_values = np.loadtxt(out_traj+'.score')
    final_unaffected_mean = np.mean(mean_values)

    number_of_times_infected = np.loadtxt(out_tg+'.intensity')
    intensity = number_of_times_infected/num_run
    for n in nodes:
        TG.add_node(n)
        TG.nodes[n]['intensity'] = intensity[n]
    TG.nodes[len(intensity)]['intensity'] = 0 #dummy node

    # load eq
    eq_dist = np.loadtxt(out_tg+'.solution')

    # plot graphs
    if CG.number_of_nodes() < PLOTLIM:
        #pass
        pos_cg = nx.kamada_kawai_layout(CG)
        pos_cg = [(pos_cg[i][0], pos_cg[i][1]) for i in range(CG.number_of_nodes())]
        pos_tg = pos_cg + [(0, 1.5 * np.max([p[1] for p in pos_cg]))]
        plot_contact_graph(CG, path.format(type='init_graph_'+type, fileformat='pdf'), init_infected, init_recovered=init_recovered, pos=pos_cg)
        plot_TG(TG, path.format(type='transm_graph_'+type, fileformat='pdf'), init_infected, eq_dist, number_of_times_infected, init_recovered=init_recovered, pos=pos_tg,  CG=CG, pos_cg=pos_cg)

    print('output done in call rust')
    return CG, TG, traj_data, final_unaffected_mean, eq_dist



def vacc_strategy_greedy(CG, init_infected, outpath, budget=2, max_steps=100, infection_rate=2.0):

    ################################################################################
    # Setup
    ################################################################################

    assert(budget>1)
    assert(outpath.endswith('.{fileformat}'))
    assert('{type}' in outpath)
    if init_infected is None:
        init_infected = [0]
    assert(len(init_infected) > 0)


    CG_reset = nx.convert_node_labels_to_integers(CG, first_label=0)  # because init is on old graph
    if set(CG.nodes()) != set(CG_reset.nodes()):
        CG = CG_reset


    ################################################################################
    # Baseline
    ################################################################################

    print('__Start vaccination baseline__')
    CG, TG, traj_data, baseline_score_novaccine, eq_dist = call_rust(CG, init_infected, path=outpath, type='baseline', infection_rate=infection_rate)

    ################################################################################
    # Performance Baselines
    ################################################################################

    # Random
    scores = list()
    candidates = [int(n) for n in CG.nodes if n not in init_infected]
    for _ in range(10):
        vaccinated = list(np.random.choice(candidates, budget))
        CG, TG, traj_data, score_rnd, eq_dist_rnd = call_rust(CG, init_infected, path=outpath, init_recovered=vaccinated,
                                                                         type='random', infection_rate=infection_rate)
        scores.append(score_rnd)
    baseline_score_random = np.mean(scores)

    # Dava
    dava_vaccinated = dava_score(CG, init_infected, budget, fast=False, inf_rate=infection_rate)
    CG, TG, traj_data, dava_baseline, eq_dist_dava = call_rust(CG, init_infected, path=outpath, type='dava', infection_rate=infection_rate, init_recovered=dava_vaccinated)

    # Dava fast
    davaf_vaccinated = dava_score(CG, init_infected, budget, fast=True, inf_rate=infection_rate)
    CG, TG, traj_dataf, davaf_baseline, eq_dist_davaf = call_rust(CG, init_infected, path=outpath, type='davafast', infection_rate=infection_rate, init_recovered=davaf_vaccinated)

    # Degree
    candidates = [n for n in CG.nodes if n not in init_infected]
    degrees = list(CG.degree(candidates))
    degrees.sort(key=lambda x:x[1], reverse=True)
    degree_vaccinated = list(map(lambda x: x[0], degrees[:budget]))
    CG, TG, traj_data, degree_baseline, eq_dist_degree = call_rust(CG, init_infected, path=outpath, type='degree', infection_rate=infection_rate, init_recovered=degree_vaccinated)

    # Between-ness
    betweenness = list(CG.betweenness_centrality())
    betweenness.sort(key=lambda x:x[1], reverse=True)
    betweenness_notinf = {node:value for (node,value) in betweenness.items() if node not in init_infected}
    betweenness_vaccinated = list(map(lambda x: x[0], betweenness_notinf[:budget]))
    CG, TG, traj_data, betweenness_baseline, eq_dist_betweenness = call_rust(CG, init_infected, path=outpath, type='betweenness', infection_rate=infection_rate, init_recovered=betweenness_vaccinated)

    # eigenvector centrality
    evcentrality = list(CG.eigenvector_centrality())
    evcentrality.sort(key=lambda x:x[1], reverse=True)
    evcentrality_notinf = {node:value for (node,value) in evcentrality.items() if node not in init_infected}
    evcentrality_vaccinated = list(map(lambda x: x[0], evcentrality_notinf[:budget]))
    CG, TG, traj_data, evcentrality_baseline, eq_dist_evcentrality = call_rust(CG, init_infected, path=outpath, type='evcentrality', infection_rate=infection_rate, init_recovered=evcentrality_vaccinated)

    # closeness centrality
    
    closeness = list(CG.closeness_centrality())
    closeness.sort(key=lambda x:x[1], reverse=True)
    closeness_notinf = {node:value for (node,value) in closeness.items() if node not in init_infected}
    closeness_vaccinated = list(map(lambda x: x[0], closeness_notinf[:budget]))
    CG, TG, traj_data, closeness_baseline, eq_dist_closeness = call_rust(CG, init_infected, path=outpath, type='closeness', infection_rate=infection_rate, init_recovered=closeness_vaccinated)


    # PageRank
    pagerank_scores = nx.pagerank(CG)
    pagerank_scores_filtered = list(filter(lambda elem: elem[0] not in init_infected, pagerank_scores.items()))
    pagerank_scores_filtered.sort(key=lambda x:x[1], reverse=True)
    pagerank_vaccinated = list(map(lambda x: x[0], pagerank_scores_filtered[:budget]))
    CG, TG, traj_data, pagerank_baseline, eq_dist_pagerank = call_rust(CG, init_infected, path=outpath, type='pagerank', infection_rate=infection_rate, init_recovered=pagerank_vaccinated)

    # Personalized PageRank
    pagerank_personalization = [1]*len(CG.nodes)
    for inf_index in [int(n) for n in CG.nodes if n in init_infected]:
        pagerank_personalization[inf_index] = 5
    pagerank_personalization_dict = dict(zip(range(0,len(pagerank_personalization)),pagerank_personalization))

    pers_pagerank_scores = nx.pagerank(CG, personalization=pagerank_personalization_dict)
    pers_pagerank_scores_filtered = list(filter(lambda elem: elem[0] not in init_infected, pers_pagerank_scores.items()))
    pers_pagerank_scores_filtered.sort(key=lambda x:x[1], reverse=True)
    pers_pagerank_vaccinated = list(map(lambda x: x[0], pers_pagerank_scores_filtered[:budget]))
    CG, TG, traj_data, pers_pagerank_baseline, eq_dist_pers_pagerank = call_rust(CG, init_infected, path=outpath, type='pers_pagerank', infection_rate=infection_rate, init_recovered=pers_pagerank_vaccinated)

    ################################################################################
    # Optimization
    ################################################################################

    results = list()
    vaccinated = list()
    used_combinations = list()
    optimization_summary = {'step_i': list(), 'score': list(), 'vaccinated': list()}
    print('__Start vaccination optimization__')
    # init
    node_candidates = [n for n in range(len(eq_dist)-1) if n not in init_infected and n not in vaccinated] #-1 important to not vacc dummy
    rank = sorted([(i,eq_dist[i]) for i in node_candidates], key= lambda  x: -x[1])
    rank = [x[0] for x in rank]
    vaccinated = rank[:1]

    first_set = None
   #step_i_helper = -1
    for step_i in range(max_steps+budget):
        step_i_str = str(step_i).zfill(6)
        removed = -1 #dummy
        print('vaccinated goes from: ', vaccinated)
        if len(vaccinated) == budget and budget > 1:
            random.shuffle(vaccinated)
            removed = vaccinated[0]
            vaccinated = vaccinated[1:]
        print('to: ', vaccinated)
        CG, TG, traj_data, score, eq_dist = call_rust(CG, init_infected, path=outpath, init_recovered=vaccinated,
                                                                         type=step_i_str, infection_rate=infection_rate)

        node_candidates = [n for n in range(len(eq_dist) - 1) if
                           n not in init_infected and n not in vaccinated and n != removed]  # -1 important to not vacc dummy
        rank = sorted([(i, eq_dist[i]) for i in node_candidates], key=lambda x: -x[1])
        #rank = [x[0] for x in rank]
        #candidate = rank[0]

        rank_best_nodes = [x[0] for x in rank[:10]]
        rank_best_p = [x[1] for x in rank[:10]]
        rank_best_normalize = np.sum(rank_best_p)
        rank_best_p = [x/rank_best_normalize for x in rank_best_p]

        for t_i in range(100): # no more than 100 tries
            if t_i == 0 or t_i == 99:
                candidate = rank[0][0]
                new_vaccinated = sorted(vaccinated + [candidate])
                if new_vaccinated not in used_combinations:
                    break
                else:
                    continue
            candidate = np.random.choice(rank_best_nodes, 1, p=rank_best_p)[0]
            new_vaccinated = sorted(vaccinated + [candidate])
            if new_vaccinated not in used_combinations:
                break

        vaccinated = sorted(new_vaccinated)

        if first_set is None and len(vaccinated) == budget:
            first_set = list(vaccinated)
        else:
            # random re-start
            if random.random() >0.95 and first_set is not None:
                vaccinated = list(first_set)

        used_combinations.append(sorted(vaccinated))
        CG, TG, traj_data, score, eq_dist = call_rust(CG, init_infected, path=outpath, init_recovered=vaccinated,
                                                                         type=step_i_str+'b', infection_rate=infection_rate)
        results.append((score, sorted(vaccinated), step_i))
        results = sorted(results, key=lambda x: -x[0])
        print('results: ', results if len(results) < 100 else results[:100])
        optimization_summary['score'].append(score)
        optimization_summary['step_i'].append(step_i)
        optimization_summary['vaccinated'].append(str(vaccinated).replace(',',';'))

        # not that doing this in each iteration might be expensive
        if step_i>2 and step_i%5==0:
            pd.DataFrame(optimization_summary).to_csv(outpath.format(type='optimization_summary', fileformat='csv'))
            plot_optimization_summary(optimization_summary, outpath.format(type='optimization_summary', fileformat='pdf'), no_vacc_score=baseline_score_novaccine, random_score = baseline_score_random, davaf=davaf_baseline, dava=dava_baseline, pagerank_baseline=pagerank_baseline, pers_pagerank_baseline=pers_pagerank_baseline,degree_baseline=degree_baseline)

    print('final results:')
    results = sorted(results, key=lambda x: -x[0])
    print(results)
    return CG, results





def analysis(CG, name, infection_rate=2.0, budget=1, max_steps=OPTIMSTEPS, init_infected=None):
    import traceback
    try:
        import os
        os.system('mkdir output/')
        os.system('mkdir output/' + name)
        CG = nx.convert_node_labels_to_integers(CG)
        outpath = 'output/'+name+'/{type}.{fileformat}'
        vacc_strategy_greedy(CG, init_infected, outpath, max_steps=max_steps, infection_rate=infection_rate, budget=budget)
    except Exception as err:
        traceback.print_tb(err.__traceback__)
        logging.info(str(err.__traceback__))
        logging.info(''.join(traceback.format_stack()))
        time.sleep(5)


# for testing (includes plotting each step)
analysis(gg.geom_graph(40), 'SmallTestNetwork', infection_rate=2.0, budget=2, init_infected=[10,11])


# Don't run experiments in test enviroment
if 'TRAVIS' not in os.environ:
    pass
    # candidates (unfort. some non-det. somewhere (maybe in graph generatino))
   # analysis(gg.geom_graph(100, 0.01), 'Exp1_Geom', infection_rate=1.0, budget=2, init_infected=[0,1,2,3,4])
   # analysis(gg.erdos_renyi(1000,0.01), 'Exp2_Erdos', infection_rate=1.3, budget=3, init_infected=[0])
    #analysis(gg.barabasi(10000), 'Exp3_BA', infection_rate=0.6, budget=10, init_infected=[100, 101,102,103,104,105,106,107,108,109])

    # speed test, might take a while
  #  speed_test()
