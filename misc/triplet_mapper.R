library(ape)
library(phytools)

args = commandArgs(trailingOnly=TRUE)


###########################
#map unbalanced triplets to species tree

e=read.tree(as.character(args[1]))
e$node.label=rep(0,length(e$tip.label)-1)
outliers=read.table('unbalanced.trp.tsv')
get_nodes_leading_to_tips<-function(tree,tips){
	tip_no=list()
	tips=as.character(tips)
	for (tip in tips){tip_no=c(tip_no,which(tree$tip.label==tip))}
	tip_no=unlist(tip_no)
	nodes=list()
	nodes=getMRCA(tree,tips)
	daughters=tree$edge[which(tree$edge[,1]==nodes[1]),2]
	#find which daughter has two target sp
	if (length(intersect(tip_no,getDescendants(tree,daughters[1]))==2)){
		#daughter1 have two ingroups
		daughter1=daughters[1]
		daughter2=daughters[2]
		single_sp=tip_no[which(tip_no %in% getDescendants(tree,daughters[2]))]
	}else{
		daughter1=daughters[2]
		daughter2=daughters[1]
	}
	#if (daughter1>length(tree$tip.label)){nodes=c(nodes,daughter1)}
	if (daughter2>length(tree$tip.label)){nodes=c(nodes,daughter2)}
	for (NOD in getDescendants(tree,daughter1)){
		if (NOD >length(tree$tip.label) & length(intersect(tip_no,getDescendants(tree,NOD)))>0){
				nodes=c(nodes,NOD)
		}
	}
	for (NOD in getDescendants(tree,nodes[1])){
		if (NOD >length(tree$tip.label) & length(intersect(single_sp,getDescendants(tree,NOD)))>0){
			nodes=c(nodes,NOD)
		}
	}
	return(nodes)
}

####################################
#count raw numbers of unbalanced triples for each node
sp_num=length(e$tip.label)
for (i in 1:length(outliers$V1)){
e$node.label[get_nodes_leading_to_tips(e,unlist(outliers[i,1:3]))-sp_num]=e$node.label[get_nodes_leading_to_tips(e,unlist(outliers[i,1:3]))-sp_num]+1
}
e$edge.length=NULL
write.tree(e,'unbalancedTriplet.raw_count.tre')

####################################
#calculating percentage of unbalanced triple for each node

trpl_totoal=e
trpl_totoal$node.label=rep(0,length(e$tip.label)-1)

for (i in (sp_num+1):(2*sp_num-1)){
	n_desced=length(getDescendants(e,i))
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
		#sample three species from ingroup
		descend_nodes=trpl_totoal$edge[which(trpl_totoal$edge[,1]==i),2]
		descend_nodes_nsp1=length(getDescendants(e,descend_nodes[1]))
		descend_nodes_nsp2=length(getDescendants(e,descend_nodes[2]))
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
write.tree(z,'unbalancedTriplet.perc.tre')