#take in one rooted species tree text with branch length in coalescent units and simulate gene trees
#output gene trees in text

geneTr_sim<-function(SpTr_text,number_geneTr){
	y=sim.coal.mpest(SpTr_text,number_geneTr)
	return(y)
}