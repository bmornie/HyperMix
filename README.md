# HyperMix

This is the official implementation accompanying

Mornie B., Colle D., Audenaert P., Pickavet M. (2026) **Flexible Generation of Random Hypergraphs Using Mixtures of Growth Mechanisms**,  
submitted to *IEEE Transactions on Network Science and Engineering*.


### Requirements
- Python 3.12+
- C++20 compatible compiler
- CMake 3.18+
- pybind11


### Installation
**1. Clone the repository**  
```bash
git clone https://github.com/bmornie/HyperMix.git
cd HyperMix
```
**2. Install Python dependencies**
```bash
pip install -r requirements.txt
```
**3. Build the C++ extension**  
Windows:
```bash
build.bat
```
Linux:
```bash
./build.sh
```

### Generating and analyzing hypergraphs
A typical workflow is as follows:  

**Optimizing weights**  
Basic usage:
```bash
python optimize.py <target_file> -l <loss_id>
```
`loss_id` can be 1 (degree EMD), 2 (degree RMSE), 3 (3-node motifs) or 4 (4-node motifs). See the paper for more info.  
  
Optional arguments:
- `-m <file>`: File of motif IDs to be included in the loss function. By default, all motifs of a given order are included. This is only used with `-l 4`.  
See `significant_motifs/` for example files. This folder contains the IDs of motifs used in the HyperMix M4<sub>5*</sub> model. The file `motif_definitions.txt` contains a representative of each motif type in a list-of-lists format.
- `-g <int>`: Maximum number of CMA-ES generations (default: 50).
- `-e <int>`: Number of candidate weights evaluated in each generation (default: 100).
- `-s <int>`: Number of hypergraphs generated for each evaluation (default: 200).
- `-f`: Toggles fixed weights. Dynamic weights are used by default.
- `-o <file>`: File in which all evaluated weights and the corresponding loss values are stored (default: `<target>_L<loss_id>.npz`).
- `-w <int>`: Number of workers used to parallellize evaluations (default: 1).  
  

**Generating Hypergraphs**  
Basic usage:
```bash
python generate.py <target_file> <weights_file> <n>
```
This will generate `n` hypergraphs using the best weights in `weights_file` (as produced by `optimize.py`).  
In principle, any npz file that stores a "data" dictionary with the "best_weights" key will do.  
  
Optional arguments:
- `-o <path>`: Output path (default: `output/`). Generated hypergraphs are stored in hyperedge list format, with filenames `<target>_<i>`, for `i` from 1 to `n`.

**Evaluating generated hypergraphs**  
Finally, we provide a script that allows to perform the evaluation from Section V.E. of the paper.  
Usage:
```bash
python analyze.py <target_file> <generated_path> <n>
```
The results are averaged over the first `n` generated hypergraphs stored in `generated_path`, assuming the naming scheme of `generate.py` is used.  
Make sure `n` is not larger than the number of stored hypergraphs!  
  
**Full example**  
The results for the Diseasome HyperMix M3 model can be reproduced with the following commands:
```bash
python optimize.py data/diseasome.txt -l 3 -g 50 -e 100 -s 200
python generate.py data/diseasome.txt diseasome_L3.npz 200
python analyze.py data/diseasome.txt output 200
```
Please note that the weights optimization will take a while with the provided settings.

### Data
The `data/` folder contains all datasets used in our experiments in the correct format.  

### Dependencies
This project uses `ankerl::unordered_dense` (included in `utilities/`), for storage of hypergraphs. We also use [our own motif counting algorithm](https://github.com/bmornie/Hypergraph-motif-counting), presented in a previous publication.  
See the corresponding repositories for details and license. No additional installation is required.

### License
This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0).  
See the `LICENSE` file for details.