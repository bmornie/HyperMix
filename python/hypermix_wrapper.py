from pathlib import Path
import sys

PYTHON_DIR = Path(__file__).resolve().parent

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))


import numpy as np
import cma
from multiprocessing.pool import ThreadPool
from functools import partial

import hypermix as hm





def generate(target, weights, max_attempts=100):
    
    if len(weights) == 10:
        dynamic_weights = True
    elif len(weights) == 5:
        dynamic_weights = False
    else:
        raise Exception("Incorrect number of weights")
    
    T = hm.read_hyperedge_list(target)
    
    Nt = T.node_count()
    Mt = T.edge_count()
    sizes = T.size_histogram()
    
    H = hm.generate_hypergraph(list(weights), Nt, Mt, sizes, dynamic_weights, max_attempts)
    
    return H


def test_weights(weights, target, loss_id, n_hypergraphs, motifs, dynamic_weights, max_attempts):
    return hm.test_weights(target, weights, loss_id, n_hypergraphs, motifs, dynamic_weights, max_attempts)

def optimize_weights(target, loss_id, max_generations, evals_per_iter, samples_per_eval, motifs, filename, dynamic_weights=True, workers=1, max_attempts=100):
    
    n_weights = 5
    
    T = hm.read_hyperedge_list(target)
    
    all_data = {"generations" : max_generations,
                "evaluations" : evals_per_iter,
                "samples" : samples_per_eval}
    
    tw = partial(test_weights, target=T, loss_id=loss_id, n_hypergraphs=samples_per_eval,
                           motifs=motifs, dynamic_weights=dynamic_weights, max_attempts=max_attempts)
    
    
    # Initial parameter set
    
    if dynamic_weights:
        initial = np.array([1]*n_weights*2)/n_weights
    else:
        initial = np.array([1]*n_weights)/n_weights
        
    es = cma.CMAEvolutionStrategy(initial, 0.5, {"popsize": evals_per_iter, "bounds": [0,1]})
    
    
    for gen in range(max_generations):
        
        weights = es.ask()
        
        # Normalize
        if dynamic_weights:
            weights = np.array([np.concatenate((w[:n_weights]/np.sum(w[:n_weights]), w[n_weights:]/np.sum(w[n_weights:]))) for w in weights])
        else:
            weights = np.array([w/np.sum(w) for w in weights])
        
        if workers > 1:
            pool = ThreadPool(workers)
            loss = pool.map(tw, weights)
        else:
            loss = np.array([tw(w) for w in weights])
        
        es.tell(weights, loss)
        
        print(
        f"Gen: {gen+1:4d} | "
        f"Evals: {es.countevals:6d} | "
        f"Best: {es.best.f:12.6g} | "
        f"Sigma: {es.sigma:.6g}",
        flush=True
        )
        
        all_data[gen+1] = {"weights" : weights,
                           "loss" : loss,
                           "best" : es.best,
                           "sigma" : es.sigma}
    
        np.savez_compressed(filename, data=all_data)
    
    
    all_data["best_weights"] = es.result.xbest
    all_data["best_loss"] = es.result.fbest
    
    np.savez_compressed(filename, data=all_data)
    
    print(f"Best: {es.result.xbest}, Loss: {es.result.fbest}")
    
    return es.result.xbest