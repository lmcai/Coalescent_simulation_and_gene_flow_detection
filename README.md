# Coalescent simulation and gene flow detection

This folder contains scripts used in Cai et al (2019) to simulate gene trees under coalescent model with bootstrap and subsequently detect gene flow using triple frequency in gene trees.
<div id="citation"></div>

<b>Citation:</b> Liming Cai, Zhenxiang Xi, Emily Moriarty Lemmon, Alan R. Lemmon, Austin Mast, Christopher E. Buddenhagen, Liang Liu, Charles C. Davis. 2019. The Perfect Storm: Gene Tree Estimation Error, Incomplete Lineage Sorting, and Ancient Gene Flow Explain the Most Recalcitrant Ancient Angiosperm Clade, Malpighiales. in review.

<b>License:</b> GPL https://www.gnu.org/licenses/gpl-3.0.html

<b>Usage contact:</b> [lcai@g.harvard.edu](mailto:lcai@g.harvard.edu)


## Dependencies

R 3.1+; R package phybase and phytools.

*If use the python script 'triple_frequency_counter.py' to count triple frequencies (suitable for species tree containing >20 terminals), python library ete3 will need to be installed.


## Input and output

<b>Input:</b> 

1. One best-estimated species tree with branch lengths measured in coalescent units (can be estimated from MPEST);

2. Bootstrap species trees with branch lengths measured in coalescent units;

3. Rooted gene trees

<b>Output:</b> 

1. tripleFr.emp.csv: Triple frequencies in the empirical gene trees;

2. BPspTr*.sim.genetrees: Bootstrapped simulated gene tree sets under multispecies coalescent model with the same amount of missing data as empirical gene trees; 

3. tripleFr.BP*.csv: Triple frequencies in the simulated gene trees. This is the expectation of triple frequency distribution under coalescent model accounting for missing data and estimation error;

4. unbalanced_triples.csv: Triples with significantly unbalanced minor frequencies;

5. unbalanced_triples_sum.tre: Species tree with node labels reflecting numbers of unbalanced triples that are associated with it.

## How to

