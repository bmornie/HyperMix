import argparse
from pathlib import Path
import numpy as np

from python.hypermix_wrapper import generate



if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("target", type=str, help="Filename of target network")
    parser.add_argument("weights", type=str, help="File of optimized weights")
    parser.add_argument("number", type=int, help="Number of hypergraphs to be generated")
    parser.add_argument("-o", "--output", type=str, default="output", metavar="PATH", help="Output path")
    
    args = parser.parse_args()
    
    
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    
    weights_file = np.load(args.weights, allow_pickle=True)
    data = weights_file["data"].item()
    weights = data["best_weights"]
    
    
    for i in range(1, args.number+1):
        H = generate(args.target, weights)
        
        with open(output_path / f"{Path(args.target).stem}_{i}.txt", 'w') as f:
            for e in H.edge_list():
                for node in e:
                    f.write(f"{node} ")
                f.write("\n")
        
        print(f"{i}/{args.number}", end="\r", flush=True)
    
    print(f"Hypergraphs stored in folder {output_path}.")