from ete3 import Tree
import itertools


#set multiple monophyletic outgroups to allow for missing data
outgroup1=["Uronema_sp", "Nephroselmis_pyriformis", "Pyramimonas_parkeae", "Monomastix_opisthostigma"]
outgroup2=["Spirogyra_sp","Chlorokybus_atmophyticus","Mesotaenium_endlicherianum"]

geneTr_text=open('').readlines()
spTr=Tree(tree_text.strip())

total_taxon=list()
for leaf in t:total_taxon.append(leaf.name)

###########################################
#initiate triple frequency dictionary
total_triple=list(itertools.permutations(taxon, 3))
triple_dict={}
#keys are the sorted species names and values are frequencies of (sp1,sp2) (sp1,sp3) (sp2,sp3) difference_of_minor_triple_frequency
for i in range(0,len(total_triple)):
		triple_dict['|'.join(sorted(total_triple[i]))]=[0,0,0,0]


def rerooter_triple_counter(i):
	print i
	t = Tree(geneTr_text[i].strip())
	taxon=list()
	for leaf in t:taxon.append(leaf.name)
	
		
	###########################################
	#gene tree rooting	
	#must have at least one outgroup
	if any(l in taxon for l in outgroup1):
		#outgroup in gene tree
		geneTr_out=[l for l in outgroup1 if l in taxon]
		if len(geneTr_out)==1:
			t.set_outgroup( t&geneTr_out[0])
		elif t.check_monophyly(values=geneTr_out, target_attr="name"):
			ancestor = t.get_common_ancestor(geneTr_out)
			t.set_outgroup(ancestor)
		else:
			#if outgroups are paraphyletic, use the monophyletic clade with largest number of outgroup species
			#raise ValueError("outgroup is paraphyletic")
			max_out=1
			print 'paraphyletic outgroup in tree' +`i`
			#annotate outgroup species
			for leaf in t:
				if leaf.name in outgroup1:
					#print leaf.name
					leaf.add_features(outgr='true')
			
			for node in t.get_monophyletic(values=["true"], target_attr="outgr"):
				if len(node)>=max_out:
					t.set_outgroup(node)
					max_out=len(node)
			
	#start to search for second outgroup clades		
	elif any(l in taxon for l in outgroup2):
		#outgroup in gene tree
		geneTr_out=[l for l in outgroup2 if l in taxon]
		if len(geneTr_out)==1:
			t.set_outgroup( t&geneTr_out[0])
		elif t.check_monophyly(values=geneTr_out, target_attr="name"):
			ancestor = t.get_common_ancestor(geneTr_out)
			t.set_outgroup(ancestor)
		else:
			#if outgroups are paraphyletic, use the monophyletic clade with largest number of outgroup species
			#raise ValueError("outgroup is paraphyletic")
			max_out=1
			print 'paraphyletic outgroup in tree' +`i`
			#annotate outgroup species
			for leaf in t:
				if leaf.name in outgroup2:
					#print leaf.name
					leaf.add_features(outgr='true')
			
			for node in t.get_monophyletic(values=["true"], target_attr="outgr"):
				if len(node)>=max_out:
					t.set_outgroup(node)
					max_out=len(node)
	
	
	else:
		print 'no outgroup in tree'+`i`
		return None


	###########################################
	#triple counting
	triples=list(itertools.permutations(taxon, 3))
	
	for j in range(0,len(triples)):
		tpl_taxa=sorted(triples[j]))
		if t.get_common_ancestor(tpl_taxa[0],tpl_taxa[2])==t.get_common_ancestor(tpl_taxa[1],tpl_taxa[2]):
			triple_dict['|'.join(sorted(total_triple[i]))][0]+=1
		elif t.get_common_ancestor(tpl_taxa[0],tpl_taxa[1])==t.get_common_ancestor(tpl_taxa[1],tpl_taxa[2]):
			triple_dict['|'.join(sorted(total_triple[i]))][1]+=1
		else:
			triple_dict['|'.join(sorted(total_triple[i]))][2]+=1
		