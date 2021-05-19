# Coalescent simulation and gene flow detection

This folder contains scripts used in Cai et al (2019) to simulate gene trees under coalescent model with bootstrap and subsequently detect gene flow using triple frequency in gene trees.
<div id="citation"></div>

<b>Citation:</b> Liming Cai, Zhenxiang Xi, Emily Moriarty Lemmon, Alan R. Lemmon, Austin Mast, Christopher E. Buddenhagen, Liang Liu, Charles C. Davis. 2019. The Perfect Storm: Gene Tree Estimation Error, Incomplete Lineage Sorting, and Ancient Gene Flow Explain the Most Recalcitrant Ancient Angiosperm Clade, Malpighiales. Systematic Biology, syaa083, https://doi.org/10.1093/sysbio/syaa083

<b>License:</b> GPL https://www.gnu.org/licenses/gpl-3.0.html

<b>Usage contact:</b> [daybreak.chua@gmail.com](mailto:daybreak.chua@gmail.com)


## Dependencies

R 3.1+; R package phybase and phytools.

Python 2 or 3; python library [ete3](http://etetoolkit.org/docs/2.3/index.html)

## Input and output

<b>Input:</b> 

1. One rooted best-estimated species tree with branch lengths measured in coalescent units (can be estimated from MPEST/ASTRAL on ML gene trees, newick format);

2. Rooted bootstrap species trees with branch lengths measured in coalescent units (can be obtained by running MPEST/ASTRAL on bootstrap gene trees, newick format);

3. Rooted gene trees, newick format.

<b>Output:</b> 

Output files will be in the working folder:

1. geneTr_sim/BP*.sim.genetrees: Simulated gene tree sets under multispecies coalescent model with the same amount of missing data as empirical gene trees;

2. [prefix].trp.csv: Triple frequencies in the empirical gene trees; 

3. geneTr_sim/BP*.trp.csv: Triple frequencies in the simulated gene trees. This is the expectation of triple frequency distribution under coalescent model accounting for missing data and estimation error;

4. unbalanced.trp.tsv: Triples with significantly unbalanced minor frequencies; 'significant' = if the differences between empirical minor frequencies is larger than the largest differences found in simulation;

5. unbalancedTriplet.sum.tre: Species tree with node labels reflecting *raw numbers* of unbalanced triples that are associated with it.

6. unbalancedTriplet.perc.tre: Species tree with node labels reflecting *percentage* of unbalanced triples that are associated with it.

## How to

Place your input files in the same folder as these scripts;

Follow three steps in `simulator_tripletCounter_tripletMapper.sh` to simulate gene trees, summarize triplet frequency distribution, and map unbalanced triplets to the species tree;

```
sh simulator_tripletCounter_tripletMapper.sh [path_to_species_tree] [path_to_bootstrap_species_trees] [path_to_gene_trees]
```

Counting triple frequencies is the most time consuming step. It can be run in parallel for bootstrap replicates with the python script `triple_frequency_counter.py`. 

## Dissecting the relative importance of gene tree estimation error, ILS, and gene flow

All of the three above mentioned factors can generate gene tree variations. Using a regression model, their relative contribution to the overall gene tree variation can be estimated.

We will use the regression method implemented in the R package [relaimpo](https://cran.r-project.org/web/packages/relaimpo/relaimpo.pdf) to assess relative importance in linear models.

The input is a matrix of dependent variable (gene tree variation) and three independent variable (gene tree error, ILS, and gene flow) across all internal nodes. Each row corresponds to these values for one internal node. The more internal nodes you have, the more powerful the regression analysis is. It is advised that there should be at least 10 observations per variable, so we ideally want at least 30 internal nodes.

<b>Input:</b> 

1. The dependent variable — Gene tree variation

It quantifies gene tree topological variation across the tree. The easiest way to calculate this value is to calculate gene concordance factor (gCF) using bootstrap gene trees in IQTREE (because it accommodates missing taxa). For exampled, if you have 20 genes, you will have 20 x 100 bootstrap gene trees assuming you conducted 100 replicates. You can combine these 2000 gene trees into one file called `BSgenetree.trees` and calculate gCF in IQTREE.

```
iqtree -t species.tre --gcf BSgenetree.trees --prefix concord
```