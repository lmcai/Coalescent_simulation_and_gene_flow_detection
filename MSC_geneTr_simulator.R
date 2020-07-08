#This is the main script that calls function to simulate gene trees from species trees, count triple frequncies, and test for gene flow based on triple frequencies. 

#load library and relevant function
library(phybase)
library(ape)
library(phytools)

###########################################################
#Input, modified the path to files
#rooted species tree with branch length in coalescent units
args = commandArgs(trailingOnly=TRUE)

# test if there are three arguments: if not, return an error
if (length(args)<3) {
  stop("Not enough arguments. Please give the path to the species tree, bootstrap species trees, and empirical gene trees.n", call.=FALSE)
}

speciesTr_path=as.character(args[1])
#rooted bootstrap species trees with branch length in coalescent units
speciesTr_BP_path=as.character(args[2])
#rooted gene trees
geneTr_path=as.character(args[3])


speciesTr=read.tree(speciesTr_path)
sp_num=length(speciesTr$tip.label)
speciesTr_BP=read.tree(speciesTr_BP_path)
print(paste('There are ',sp_num,' species.',sep=''))
geneTr=read.tree(geneTr_path)
geneTr_num=length(geneTr)
print(paste('There are ',geneTr_num,' gene trees.',sep=''))

###########################################################
#simulate gene trees under coalescent model for BP species trees

geneTr_sim<-function(SpTr_text,number_geneTr){
	y=sim.coal.mpest(SpTr_text,number_geneTr)
	return(y)
}

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
	write(paste(geneTr_sim_text,collapse='\n'),paste('./geneTr_sim/',i,'.tem.genetrees',sep=''))
	simulated_gene_trees=read.tree(paste('./geneTr_sim/',i,'.tem.genetrees',sep=''))
	
	#simulate missing data according to empirical data
	for (sp in speciesTr$tip.label){
		num_missing=geneTr_num-as.numeric(sp_missing[sp])
		tr2prune=sample(1:geneTr_num,num_missing)
		for (j in tr2prune){
			simulated_gene_trees[[j]]=drop.tip(simulated_gene_trees[[j]],sp)
		}
	}
	
	
	#output to file
	write.tree(simulated_gene_trees,paste('./geneTr_sim/BP',i,'.sim.genetrees',sep=''))
	
}
