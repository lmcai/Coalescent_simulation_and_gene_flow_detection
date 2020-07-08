from ete3 import Tree
import itertools
import sys

#set multiple monophyletic outgroups to allow for missing data

geneTr_text=open(sys.argv[1]).readlines()
spTr_text=open(sys.argv[2]).readline()

spTr = Tree(spTr_text.strip())
total_taxon=list()
for leaf in spTr:total_taxon.append(leaf.name)

###########################################
#initiate triple frequency dictionary
total_triple=list(itertools.combinations(total_taxon, 3))
triple_dict={}
#keys are the sorted species names and values are frequencies of (sp1,sp2) (sp1,sp3) (sp2,sp3) difference_of_minor_triple_frequency
for i in range(0,len(total_triple)):
		triple_dict['|'.join(sorted(total_triple[i]))]=[0,0,0,0]

###########################################
#loop though genes 
for i in range(0,len(geneTr_text)):
	print 'Processing tree'+`i`
	t = Tree(geneTr_text[i].strip())
	taxon=list()
	for leaf in t:taxon.append(leaf.name)
	###########################################
	#gene tree rooting	
	#ancestor=t.get_common_ancestor('Clus','Cdub')
	try:
	#	t.set_outgroup(ancestor)
	###########################################
	#triple counting
		triples=list(itertools.combinations(taxon, 3))
	
		for j in range(0,len(triples)):
			tpl_taxa=sorted(triples[j])
			if t.get_common_ancestor(tpl_taxa[0],tpl_taxa[2])==t.get_common_ancestor(tpl_taxa[1],tpl_taxa[2]):
				triple_dict['|'.join(tpl_taxa)][0]+=1
			elif t.get_common_ancestor(tpl_taxa[0],tpl_taxa[1])==t.get_common_ancestor(tpl_taxa[1],tpl_taxa[2]):
				triple_dict['|'.join(tpl_taxa)][1]+=1
			else:
				triple_dict['|'.join(tpl_taxa)][2]+=1
	except:
		print('bad rooting')
		pass

out=open(sys.argv[1].split('.')[0]+'.trp.tsv','a')

for k, v in triple_dict.iteritems():
	out.write(k+'\t'+'\t'.join([str(l) for l in v])+'\n')

out.close()

