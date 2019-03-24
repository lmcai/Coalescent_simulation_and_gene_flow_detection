#This R script is used to count the triple frequencies in a collection of gene trees
#The required input is one species tree (newick form in one file) and all of your rooted gene trees (newick form in one file)
#Liming Cai (lcai@g.harvard.edu)



count_triple_treq<-function(phy,matrix){
	for (i in 1:length(matrix[1,])){
		if (all(matrix[1:3,i] %in% phy$tip.label)){
			pruned.tree<-drop.tip(phy,phy$tip.label[-match(matrix[1:3,i], phy$tip.label)])
			if (is.monophyletic(pruned.tree,matrix[c(1,2),i])){
				matrix[4,i]=as.numeric(matrix[4,i])+1
				next	
			}
			if (is.monophyletic(pruned.tree,matrix[c(1,3),i])){
				matrix[5,i]=as.numeric(matrix[5,i])+1
				next
			}
			if (is.monophyletic(pruned.tree,matrix[c(2,3),i])){
				matrix[6,i]=as.numeric(matrix[6,i])+1
				next
			}
		}
	}
	return(matrix)
}



triple_writer<-function(speciesTr,MLgeneTr,output_nam){
	####read in a species tree and initiate values in triple frequency
	spNam=speciesTr$tip.label
	ML_triples=combn(spNam,3)
	ML_triples=data.frame(ML_triples,stringsAsFactors=FALSE)
	ML_triples=rbind(ML_triples,rep(0,length(ML_triples[1,])))
	ML_triples=rbind(ML_triples,rep(0,length(ML_triples[1,])))
	ML_triples=rbind(ML_triples,rep(0,length(ML_triples[1,])))

	##########calculate triple frequency for gene trees########
	for (j in 1:length(MLgeneTr)){
		print(paste('Processing tree',j))
		ML_triples=count_triple_treq(MLgeneTr[[j]],ML_triples)
	}
	write.csv(ML_triples,output_nam)
}


