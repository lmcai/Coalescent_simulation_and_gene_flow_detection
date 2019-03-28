#This is the main script that calls function to simulate gene trees from species trees, count triple frequncies, and test for gene flow based on triple frequencies. 

#load library and relevant function
library(phybase)
library(ape)
library(phytools)

source('./count_triple_refq.R')
source('./geneTr_simulator.R')


#input
#species tree with branch length in coalescent units
speciesTr_path=
speciesTr=read.tree(speciesTr_path)
sp_num=length(speciesTr$tip.label)

#bootstrap species tree with branch length in coalescent units
speciesTr_BP_path=
speciesTr_BP=read.tree(speciesTr_BP_path)

#gene trees
geneTr_path=
geneTr=read.tree(geneTr_path)
geneTr_num=length(geneTr)

#output
out_dir=

###########################################################
#count triple frequency in empirical data
triple_writer(speciesTr,geneTr,paste(out_dir,'/tripleFr.emp.csv',sep=''))

#output triple frequency matrix written to 'out_dir/tripleFr.emp.csv'

###########################################################
#simulate gene trees under coalescent model for BP species trees

#prepare list of missing data for each species
sp_missing=list()
#initiate
for (sp in speciesTr$tip.label){sp_missing[sp]=0}
for (sp in speciesTr$tip.label){
	for (tr in geneTr){
		if (sp %in% tr$tip.label){sp_missing[sp]=as.numeric(sp_missing[sp])+1}
	}
}


speciesTr_BP_text=readLines(speciesTr_BP_path)
for (i in 1:length(speciesTr_BP_text)){
	geneTr_sim_text=geneTr_sim(speciesTr_BP_text[i],geneTr_num)
	simulated_gene_trees=read.tree(geneTr_sim_text)
	
	#simulate missing data according to empirical data
	for (sp in speciesTr$tip.label){
		num_missing=geneTr_num-as.numeric(sp_missing[sp])
		tr2prune=sample(1:geneTr_num,num_missing)
		for (j in tr2prune){
			simulated_gene_trees[[j]]=drop.tip(simulated_gene_trees[[j]],sp)
		}
	}
	
	
	#output to file
	write.tree(simulated_gene_trees,paste(out_dir,'/geneTr_sim/BPspTr',i,'.sim.genetrees',sep=''))
	
	#count triple frequency in the simulated data
	triple_writer(speciesTr,simulated_gene_trees,paste(out_dir,'/tripleFr_sim/tripleFr.BP',i,'.csv',sep=''))
}

###########################################################
#compare the empirical and the simulated triple frequency distribution and identify significantly unbalanced triples

tripleFr_outlier_detect

outliers





#output to csv
write.csv(outliers,paste(out_dir,'/unbalanced_triples.csv',sep=''))

###########################################################
#summarize the percentage of unbalanced triples onto species trees

#read unrooted mpest species tree
e=read.tree('/Users/lcai/Downloads/Malpighiales_ILS/coal_ana/ExonIntron/total423.mpest.mb.95cs.sum.tre')
e$node.label=rep(0,length(e$node.label))

get_nodes_leading_to_tips<-function(tree,tips){
	tip_no=list()
	for (tip in tips){tip_no=c(tip_no,which(tree$tip.label==tip))}
	tip_no=unlist(tip_no)
	nodes=list()
	nodes=getMRCA(tree,tips)
	for (NOD in getDescendants(tree,nodes[1])){
		if (NOD >length(tree$tip.label) & length(intersect(tip_no,getDescendants(tree,NOD)))>0){
			nodes=c(nodes,NOD)
		}
	}
	return(nodes)
}

#count raw numbers of unbalanced triples for each node
for (i in 2:length(outliers[1,])){
	e$node.label[get_nodes_leading_to_tips(e,as.character(outliers[1:3,i]))-sp_num]=e$node.label[get_nodes_leading_to_tips(e,as.character(outliers[1:3,i]))-sp_num]+1
}
e$edge.length=NULL
write.tree(e,paste(out_dir,'/unbalanced_triples_sum.tre',sep=''))

