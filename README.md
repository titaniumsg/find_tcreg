# T-CREGs

#### For a detailed description of the methods, please refer to the [article](https://doi.org/10.1101/2024.06.04.597418):
Gupta S, Sgourakis NG. A structure-guided approach to predict MHC-I restriction of T cell receptors for public antigens. (2024). bioRxiv. doi.org/10.1101/2024.06.04.597418.

Corresponding author: Nikolaos G. Sgourakis, nikolaos.sgourakis@pennmedicine.upenn.edu

### System requirements:
    - MacOS/Linux (tested on MacOS 14.5)
    - Python 3.8+ (tested on 3.8.15)
    - Anaconda (tested on 4.11.0)
        - biopython
        - pymol

### Installation:
1. Create a conda environment with Python 3.8+, [biopython](https://anaconda.org/conda-forge/biopython), and [PyMOL](https://pymol.org/conda/) installed.
2. Download the repository.

### Usage

The code has two modes. The first mode assumes the given receptor's binding mode is similar to that of TCRs and utilizes the previously identified TCR-contacting residues from TCR:pHLA structures in TCR3d. The second mode recalculates the T-CREGs using the receptor:pHLA complex structure as input, taking into account the binding mode in the identification of contact residues. 

The default T-CREG file is the "common" 21 T-CREGs. If you'd like to use all T-CREGs, include the following flag: `-tcreg_file tcregs_greedy_no_cutoff.txt`. Custom T-CREG files can also be inputted using this flag.

#### Mode 1

This mode only requires the set of HLAs anticipated to bind the peptide and the output directory. An example of the set of HLAs can be found in the `example` directory and can be run using the following command:

`python tcreg_pipeline.py -hla_file example/PHOX2B_strong_binding_alleles.txt -output_dir example/`

#### Mode 2

This mode requires everything Mode 1 requires as well as the path to the structure (`-structure`), the HLA chain name (`-hla_chain`), and the receptor chain(s) (`-receptor_chain`). If there are multiple receptor chains, they can be separated by spaces e.g., `-receptor_chain D E F G`. The example of mode 2 can be run using the following command:

`python tcreg_pipeline.py -hla_file example_with_structure/PHOX2B_strong_binding_alleles.txt -output_dir example_with_structure/ -structure example_with_structure/10LH.pdb -hla_chain A -receptor_chain H L`

