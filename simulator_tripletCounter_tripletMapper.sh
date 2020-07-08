#0. input files, please see read me for the format of these files
#path to rooted species tree (newick format), branch length in coalescent units
speciesTr=$1
#path to rooted bootstrap species trees (newick format), branch length in coalescent units
speciesTr_BP=$2
#path to rooted empirical gene trees (newick format)
geneTr=$3

#1. Simulate gene trees under coalescent model for each BP species tree
#variations in the resulting simulated gene trees will reflect estimation error, ILS, and missing data in the empirical data set
printf "Simulating gene trees under MSC model...\n"
mkdir geneTr_sim
Rscript --vanilla MSC_geneTr_simulator.R $speciesTr $speciesTr_BP $geneTr
rm geneTr_sim/*.tem.genetrees
#excepted result: sets of simulated gene trees in the folder geneTr_sim

#2. Count triplet frequency in empirical gene trees and simulated gene trees
# This can take hours if there are >30 species, so it is strongly advised to distribute the work to the cluster. For example, submit a job to count triplet frequency for one set of gene trees.
printf "Counting triplet frequency in the empirical data...\n"
python triple_frequency_counter.py $geneTr $speciesTr
#excepted result: *.trp.tsv
#format: column1--species names of the triplet, sorted alphabetically; column2--triplet frequencies of (sp1,sp2);column3--triplet frequencies of (sp1,sp3);column4--triplet frequencies of (sp2,sp3).

printf "Now you decide how to count triplet frequency for the simulated bootstrap replicates. !!This script will skip this step!! You can distribute them into separate jobs on cluster...\n"
printf "example:\npython triple_frequency_counter.py geneTr_sim/BP1.sim.genetrees $speciesTr\npython triple_frequency_counter.py geneTr_sim/BP2.sim.genetrees $speciesTr\npython triple_frequency_counter.py geneTr_sim/BP3.sim.genetrees $speciesTr\n......\n"

#3. Find significantly unbalanced triplets and map to species tree
python find_unbalanced_triplets.py $speciesTr
Rscript --vanilla triplet_mapper.R $speciesTr