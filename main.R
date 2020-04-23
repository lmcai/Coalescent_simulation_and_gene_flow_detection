#This is the main script that calls function to simulate gene trees from species trees, count triple frequncies, and test for gene flow based on triple frequencies. 

#load library and relevant function
library(phybase)
library(ape)
library(phytools)

source('./count_triple_refq.R')
source('./geneTr_simulator.R')


#input
#species tree with branch length in coalescent units
speciesTr_path='./FILE'
speciesTr=read.tree(speciesTr_path)
sp_num=length(speciesTr$tip.label)

#bootstrap species trees with branch length in coalescent units
speciesTr_BP_path='./FILE'
speciesTr_BP=read.tree(speciesTr_BP_path)

#gene trees
geneTr_path='./FILE'
geneTr=read.tree(geneTr_path)
geneTr_num=length(geneTr)

#output
out_dir='./FOLDER'

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
	#simulate gene trees under coalescent model
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
#read in csvs of triples frequencies form the empirical data as well as the simulated data
#output triples with differences in minor triple frequency exceed the maximum value seen in simulation

#output to csv
write.csv(outliers,paste(out_dir,'/unbalanced_triples.csv',sep=''))

###########################################################
#summarize the percentage of unbalanced triples onto species trees

#read unrooted mpest species tree
e=speciesTr
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
write.tree(e,paste(out_dir,'/unbalanced_triples_raw_count.tre',sep=''))

############
#calculating percentage of unbalanced triple for each node
#######
trpl_totoal=e
trpl_totoal$node.label=rep(0,length(e$tip.label)-1)

for (i in (sp_num+1):(2*sp_num-1)){
	n_desced=length(getDescendants(x,i))
	sister_desced=sp_num-n_desced
	if (sister_desced>0){
		#sample 2 ingroups and one outgroup
		trpl_totoal$node.label[i-sp_num]=trpl_totoal$node.label[i-sp_num]+n_desced*(n_desced-1)/2*sister_desced
	}
	if (sister_desced>1){
		#sample one ingroup and two outgroups
		trpl_totoal$node.label[i-sp_num]=trpl_totoal$node.label[i-sp_num]+n_desced*sister_desced*(sister_desced-1)/2
	}
	if (n_desced>2){
		#sample tree species from ingroup
		descend_nodes=trpl_totoal$edge[which(trpl_totoal$edge[,1]==i),2]
		descend_nodes_nsp1=length(getDescendants(x,descend_nodes[1]))
		descend_nodes_nsp2=length(getDescendants(x,descend_nodes[2]))
		if (descend_nodes_nsp1>1){
			trpl_totoal$node.label[i-sp_num]=trpl_totoal$node.label[i-sp_num]+descend_nodes_nsp1*(descend_nodes_nsp1-1)/2*descend_nodes_nsp2
		}
		if (descend_nodes_nsp2>1){
			trpl_totoal$node.label[i-sp_num]=trpl_totoal$node.label[i-sp_num]+descend_nodes_nsp2*(descend_nodes_nsp2-1)/2*descend_nodes_nsp1
		}
	}
}
z=e
z$node.label=as.integer(e$node.label)/as.integer(trpl_totoal$node.label)
write.tree(z,paste(out_dir,'/unbalanced_triples_perc.tre',sep=''))