#This is the main script that calls function to simulate gene trees from species trees, count triple frequncies, and test for gene flow based on triple frequencies. 

#load library and relevant function
library(phybase)
library(ape)
library(phytools)

source('./count_triple_refq.R')
source('./geneTr_simulator.R')


#input

speciesTr_path=
speciesTr=read.tree(speciesTr_path)

speciesTr_BP_path=
speciesTr_BP=read.tree(speciesTr_BP_path)

geneTr_path=
geneTr=read.tree(geneTr_path)
geneTr_num=length(geneTr)

#output
out_prefix=

###########################################################
#count triple frequency in empirical data
triple_writer(speciesTr,geneTr,paste(out_prefix,'_tripleFr.emp.csv',sep=''))

#output triple frequency matrix written to 'out_prefix__triple_fr_emp.csv'

###########################################################
#simulate gene trees under coalescent model for BP species trees

#prepare dictionary of missing data for each species





speciesTr_BP_text=readLines(speciesTr_BP_path)
for (i in 1:length(speciesTr_BP_text)){
	geneTr_sim_text=geneTr_sim(speciesTr_BP_text[i],geneTr_num)
	simulated_gene_trees=read.tree(geneTr_sim_text)
	
	#simulate missing data according to empirical data
	
	
	
	
	
	
	
	
	
	#output to file
	write.tree(simulated_gene_trees,paste(out_prefix,'.BPspTr',i,'.sim.genetrees',sep=''))
	
	#count triple frequency in the simulated data
	triple_writer(speciesTr,simulated_gene_trees,paste(out_prefix,'_tripleFr.BP',i,'.csv',sep=''))
}

###########################################################
#compare the empirical and the simulated triple frequency distribution and identify significantly unbalanced triples


#read triple frequency and find significant unbalanced minor triple through chisq test pvalue
#x=read.csv('total423.MLgene.triFr.csv')
#chisq test
#y=list()
#for (i in 2:length(x[1,])){
#	ordered_num=sort(as.numeric(as.character(x[4:6,i])))
#	y=c(y,chisq.test(ordered_num,p=c((1-ordered_num[3]/sum(ordered_num))/2,(1-ordered_num[3]/sum(ordered_num))/2,ordered_num[3]/sum(ordered_num)))$p.value)
#}
#a=which(y<0.05)+1















#summarize the percentage of unbalanced triples onto species trees

a=read.csv('/Users/lcai/Downloads/Malpighiales_ILS/Triple_frequency/unbalance_minor_triples.csv')
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

#
for (i in 2:length(a[1,])){
	e$node.label[get_nodes_leading_to_tips(e,as.character(a[1:3,i]))-65]=e$node.label[get_nodes_leading_to_tips(e,as.character(a[1:3,i]))-65]+1
}

e$edge.length=NULL
write.tree(e,'/Users/lcai/Downloads/Malpighiales_ILS/Triple_frequency/Speciestree_unbalanced_triples.tre')

