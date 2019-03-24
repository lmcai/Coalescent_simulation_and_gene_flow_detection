x=read.csv('/Users/lcai/Downloads/Malpighiales_ILS/Triple_frequency/total423.MLgene.triFr.csv',header=F)

outlier=data.frame(x[,2],stringsAsFactors = F)

for (i in 2:length(x[1,])){
	#print(i)
	ordered_num=sort(as.numeric(as.character(x[5:7,i])))
	m1=ordered_num[1]
	m2=ordered_num[2]
	
	#empirical observation to test
	test <- matrix(c(m1,m2-m1),byrow = TRUE,ncol = 2, nrow = 1)
	if (test[1]<10){
		next
	}else{
	#read in simulations adn its 95% ellipse
	a=read.csv(paste(x[1,i],'.txt',sep=''),header=F)
	p=list()
	for (j in 1:length(a$V1)){
		p=c(p,min(as.numeric(c(a$V5[j],a$V6[j]))),a$V7[j])
	}
	p=as.integer(p)
	p=matrix(p,byrow = TRUE,ncol = 2)
	
	if (test[1]<quantile(p[,2],0.1) | test[1]>quantile(p[,2],0.9)){
		if (test[2]>max(a$V7)){
		#if (test[2]>quantile(p[,2],0.95)){
			outlier=cbind(outlier,x[,i])
		}
	}else{
		Z <- pointsToEllipsoid(test, cov(p), colMeans(p))
	
		inside <- ellipseInOut(Z, p = 0.95)
		if (length(inside==TRUE)==0){
			outlier=cbind(outlier,x[,i])
		}
	}
	}
}

write.csv(outlier,'total423.unbalance_minor_triples.csv')
