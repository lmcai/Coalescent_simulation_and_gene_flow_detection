###########################################################
#compare the empirical and the simulated triple frequency distribution and identify significantly unbalanced triples
#read in csvs of triples frequencies form the empirical data as well as the simulated data
#output triples with differences in minor triple frequency exceed the maximum value seen in simulation

#output to csv
write.csv(outliers,paste(out_dir,'/unbalanced_triples.csv',sep=''))

###########################################################
#summarize the number of unbalanced triples onto species trees

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
	n_desced=length(getDescendants(trpl_totoal,i))
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
		descend_nodes_nsp1=length(getDescendants(trpl_totoal,descend_nodes[1]))
		descend_nodes_nsp2=length(getDescendants(trpl_totoal,descend_nodes[2]))
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