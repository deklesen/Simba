import random
import networkx as nx
from evaluation import vacc_strategy_greedy
import logging

def baseline_Simba(CG, init_infected, budget, outpath, infection_rate, max_steps, **kwargs):
    import traceback
    try:
        
        return vacc_strategy_greedy(CG, init_infected, outpath, max_steps=max_steps, infection_rate=infection_rate, budget=budget, Simple=False)
    except Exception as err:
        traceback.print_tb(err.__traceback__)
        logging.info(str(err.__traceback__))
        logging.info(''.join(traceback.format_stack()))
        exit()
       # time.sleep(5)

def baseline_SimbaSimple(CG, init_infected, budget, outpath, infection_rate, max_steps, **kwargs):
    import traceback
    try:
        
        return vacc_strategy_greedy(CG, init_infected, outpath, max_steps=max_steps, infection_rate=infection_rate, budget=budget, Simple=True)
    except Exception as err:
        traceback.print_tb(err.__traceback__)
        logging.info(str(err.__traceback__))
        logging.info(''.join(traceback.format_stack()))
        exit()
       # time.sleep(5)



def baseline_none(*args,**kwargs):
    return []

def baseline_random(CG, init_infected, budget, **kwargs):
    candidates = [int(n) for n in CG.nodes if n not in init_infected]
    return random.choices(candidates, k = budget)


def dava_score(CG, init_infected, budget, fast=False, inf_rate=None):
    from standalone_dava import dava_intervention
    edge_list = list(CG.edges())
    if inf_rate is None:
        si_transmit_prob = 1.0
    else:
        si_transmit_prob = inf_rate/(inf_rate + 1)
    edge_list = [(x,y,0.5) for x,y in edge_list]
    node_ids = dava_intervention(edge_list, infected_list=init_infected, recovered_list=[], k=budget, plotting=False, fast=fast)
    if not (len(node_ids) == budget):
        print(len(node_ids) , budget)
    assert(len(node_ids) <= budget)
    for node in node_ids:
        assert(node not in init_infected)
    return node_ids

def baseline_dava(CG, init_infected, budget, **kwargs):
    dava_vaccinated = dava_score(CG, init_infected, budget, fast=False)
    return dava_vaccinated

def baseline_davaf(CG, init_infected, budget,**kwargs):
    davaf_vaccinated = dava_score(CG, init_infected, budget, fast=True)
    return davaf_vaccinated

def baseline_degree(CG, init_infected, budget,**kwargs):
    # Degree
    candidates = [n for n in CG.nodes if n not in init_infected]
    degrees = list(CG.degree(candidates))
    degrees.sort(key=lambda x:x[1], reverse=True)
    degree_vaccinated = list(map(lambda x: x[0], degrees[:budget]))
    return degree_vaccinated

def baseline_betweenness(CG, init_infected, budget,**kwargs):

    # Between-ness
    betweenness = sorted(nx.betweenness_centrality(CG).items(), key=lambda x:x[1], reverse=True)
    
    #betweenness.sort(key=lambda x:x[1], reverse=True)
    betweenness_notinf = [(node,value) for (node,value) in betweenness if node not in init_infected]
    betweenness_vaccinated = list(map(lambda x: x[0], betweenness_notinf[:budget]))
    return betweenness_vaccinated
    

def baseline_EVcentrality(CG, init_infected, budget,**kwargs):
# Between-ness
    evc = sorted(nx.eigenvector_centrality(CG, max_iter=500, tol=1e-4).items(), key=lambda x:x[1], reverse=True)
    
    #evc.sort(key=lambda x:x[1], reverse=True)
    evc_notinf = [(node,value) for (node,value) in evc if node not in init_infected]
    evc_vaccinated = list(map(lambda x: x[0], evc_notinf[:budget]))
    return evc_vaccinated

def baseline_closenessCentrality(CG, init_infected, budget,**kwargs):
    close = sorted(nx.closeness_centrality(CG).items(), key=lambda x:x[1], reverse=True)
    
    #close.sort(key=lambda x:x[1], reverse=True)
    close_notinf = [(node,value) for (node,value) in close if node not in init_infected]
    close_vaccinated = list(map(lambda x: x[0], close_notinf[:budget]))
    return close_vaccinated


def baseline_Pagerank(CG, init_infected, budget,**kwargs):
    pagerank_scores = nx.pagerank(CG)
    pagerank_scores_filtered = list(filter(lambda elem: elem[0] not in init_infected, pagerank_scores.items()))
    pagerank_scores_filtered.sort(key=lambda x:x[1], reverse=True)
    pagerank_vaccinated = list(map(lambda x: x[0], pagerank_scores_filtered[:budget]))
    return pagerank_vaccinated

def baseline_PersPagerank(CG, init_infected, budget,**kwargs):
    # Personalized PageRank
    pagerank_personalization = [1]*len(CG.nodes)
    for inf_index in [int(n) for n in CG.nodes if n in init_infected]:
        pagerank_personalization[inf_index] = 10
    pagerank_personalization_dict = dict(zip(range(0,len(pagerank_personalization)),pagerank_personalization))

    pers_pagerank_scores = nx.pagerank(CG, personalization=pagerank_personalization_dict)
    pers_pagerank_scores_filtered = list(filter(lambda elem: elem[0] not in init_infected, pers_pagerank_scores.items()))
    pers_pagerank_scores_filtered.sort(key=lambda x:x[1], reverse=True)
    pers_pagerank_vaccinated = list(map(lambda x: x[0], pers_pagerank_scores_filtered[:budget]))
   
    return pers_pagerank_vaccinated