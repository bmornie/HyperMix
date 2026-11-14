import argparse
from pathlib import Path
import numpy as np

import python.hypermix as hm



def relative_deviation(real, generated):
    distances = []
    for g in generated:
        distances.append(np.abs(real - g) / np.abs(real))
    return np.mean(distances)

def motif_MRE(real, generated): # Mean relative error
    distances = []
    for g in generated:
        MARE = 1/len(real) * np.sum(np.abs(np.array(real) - np.array(g)) / (np.array(real)+1))
        distances.append(MARE)
    return np.mean(distances)

def wasserstein1(real, generated):
    h1 = np.asarray(real, dtype=float)
    
    distances = []
    for g in generated:
        h2 = np.asarray(g, dtype=float)
    
        # Normalize to probability distributions
        p = h1 / h1.sum()
        q = h2 / h2.sum()
    
        # Pad to the same length (missing degrees have zero count)
        n = max(len(p), len(q))
        p = np.pad(p, (0, n - len(p)))
        q = np.pad(q, (0, n - len(q)))
        
        # W1 for integer-valued degrees with unit spacing
        distances.append(float(np.abs(np.cumsum(p - q)[:-1]).sum()))
    
    return np.mean(distances)



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("target", type=str, help="Filename of target network")
    parser.add_argument("generated", type=str, help="Path to generated hypergraph files")
    parser.add_argument("number", type=int, help="Number of generated hypergraphs to use in computation of average deviations")
    
    args = parser.parse_args()
    
    
    T = hm.read_hyperedge_list(args.target)
    
    
    print("Computing properties of generated hypergraphs...", flush=True)
    
    degree = []
    intersect = []
    M3 = []
    diameter = []
    for i in range(1, args.number+1):
        H = hm.read_hyperedge_list(f"{args.generated}/{Path(args.target).stem}_{i}.txt")
        degree.append(H.degree_histogram())
        M3.append(H.motif_count(3))
        intersect.append(H.intersection_histogram())
        diameter.append(hm.effective_diameter(H))
        
        print(f"{i}/{args.number}", end="\r", flush=True)
    
    
    dev_degree = wasserstein1(T.degree_histogram(), degree)
    dev_intersect = wasserstein1(T.intersection_histogram(), intersect)
    dev_M3 = motif_MRE(T.motif_count(3), M3)
    dev_diameter = relative_deviation(hm.effective_diameter(T), diameter)
    
    
    print(f"Degree: {dev_degree:.3f}")
    print(f"Intersection size: {dev_intersect:.3f}")
    print(f"3-node motifs: {dev_M3:.3f}")
    print(f"Effective diameter: {dev_diameter:.3f}")