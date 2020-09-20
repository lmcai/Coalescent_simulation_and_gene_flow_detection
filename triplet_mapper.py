from ete3 import Tree
import sys
#speciesTr=Tree('ML.mpest.rooted.tre',format=1)

speciesTr=Tree(sys.argv[1],format=1)
for nd in speciesTr.traverse():
	if not nd.is_leaf():nd.name=0
#outliers=open('outliers.filtered.tsv').readlines()
#triple_frequency=open('ML.trp.tsv').readlines()
outliers=open(sys.argv[2]).readlines()
triple_frequency=open(sys.argv[1].split('.')[0]+'.trp.tsv').readlines()

#dictionary of the second most frequent topology in a triplet, candidates of introgression
intro_branch_in_triplet={}
for i in triple_frequency:
	a=i.split()
	b=[int(j) for j in a[1:4]]
	b.sort()
	if int(a[1])==b[1]:
		#(sp1,sp2) introgression 
		intro_branch_in_triplet[a[0]]=a[0].split('|')[:2]
	elif int(a[2])==b[1]:
		#(sp1,sp3) introgression 
		intro_branch_in_triplet[a[0]]=[a[0].split('|')[0],a[0].split('|')[2]]
	elif int(a[3])==b[1]:
		#(sp2,sp3) introgression 
		intro_branch_in_triplet[a[0]]=a[0].split('|')[1:]
		

#only get nodes that are involved in introgression
def get_introgression_nodes_leading_to_tips(tr,tips,introgression_tips):
	mrca_all=tr.get_common_ancestor(tips)
	mrca_introgression=tr.get_common_ancestor(introgression_tips)
	if not mrca_all==mrca_introgression:
		#introgression take place in the two closely related species
		return(tr)
	else:
		#add counts to all nodes in introgression branches
		non_introgression_tip=[j for j in tips if not j in introgression_tips]
		#((sp,sp)node1,sp)node2
		if tr.get_common_ancestor(non_introgression_tip[0],introgression_tips[0])!=mrca_all:
			#node1
			node1=tr.get_common_ancestor(non_introgression_tip[0],introgression_tips[0])
			node1.name=node1.name+1
			for nd in node1.get_descendants():
				if introgression_tips[0] in [leaf.name for leaf in nd] and not nd.is_leaf():
					nd.name=nd.name+1
			#node2
			node2=mrca_all
			node2.name=node2.name+1
			for nd in node2.get_descendants():
				if introgression_tips[1] in [leaf.name for leaf in nd] and not nd.is_leaf():
					nd.name=nd.name+1
			return(tr)
		elif tr.get_common_ancestor(non_introgression_tip[0],introgression_tips[1])!=mrca_all:
			#node1
			node1=tr.get_common_ancestor(non_introgression_tip[0],introgression_tips[1])
			node1.name=node1.name+1
			for nd in node1.get_descendants():
				if introgression_tips[1] in [leaf.name for leaf in nd] and not nd.is_leaf():
					nd.name=nd.name+1
			#node2
			node2=mrca_all
			node2.name=node2.name+1
			for nd in node2.get_descendants():
				if introgression_tips[0] in [leaf.name for leaf in nd] and not nd.is_leaf():
					nd.name=nd.name+1
			return(tr)

for l in outliers:
	speciesTr=get_introgression_nodes_leading_to_tips(speciesTr,l.split()[:3],intro_branch_in_triplet['|'.join(l.split()[:3])])

speciesTr.write(outfile='unbalanced_triples_raw_count.tre',format=1)

#####################
#calculate Reticulation In dex for each node (raw count/total number of triplets for each node)



	