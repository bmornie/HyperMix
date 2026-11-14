import argparse
from pathlib import Path

from python.hypermix_wrapper import optimize_weights



if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("target", type=str, help="Filename of target network")
    parser.add_argument("-l", "--loss", type=int, required=True, metavar="N", help="Id of loss function (1: degree EMD, 2: degree RMSE, 3: 3-node motifs, 4: 4-node motifs)")
    parser.add_argument("-m", "--motifs", type=str, metavar="FILE", help="File containing ids of 4-node motifs. Only used with --loss 4")
    parser.add_argument("-g", "--gen", type=int, default=50, metavar="N", help="Number of generations")
    parser.add_argument("-e", "--eval", type=int, default=100, metavar="N", help="Evaluations per generation")
    parser.add_argument("-s", "--samples", type=int, default=200, metavar="N", help="Samples per evaluation")
    parser.add_argument("-f", "--fixed", action="store_true", help="Use fixed instead of dynamic weights")
    parser.add_argument("-o", "--output", type=str, metavar="FILE", help="Output filename")
    parser.add_argument("-w", "--workers", type=int, default=1, metavar="N", help="Workers used for parallel evaluations")
    
    args = parser.parse_args()
    
    if args.loss not in (1,2,3,4):
        raise Exception(f"Loss id {args.loss} is invalid. Please use 1, 2, 3 or 4.")
    
    if args.output is not None:
        output = args.output + ".npz"
    else:
        output = f"{Path(args.target).stem}_L{args.loss}.npz"
    
    
    if args.motifs is None:
        motifs = []
    else:
        with open(args.motifs, "r") as file:
            motifs = [int(line.strip()) for line in file]
    
    
    
    
    optimize_weights(args.target,
                     args.loss,
                     args.gen,
                     args.eval,
                     args.samples,
                     motifs,
                     output,
                     dynamic_weights = not args.fixed,
                     workers = args.workers
                     )