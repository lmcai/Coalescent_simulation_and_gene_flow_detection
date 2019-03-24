#This is the main script that calls function to simulate gene trees from species trees, count triple frequncies, and test for gene flow based on triple frequencies. 

#load library and relevant function
library(phybase)
library(ape)

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
triple_writer(speciesTr,geneTr,out_prefix)

#output triple frequency matrix written to 'out_prefix__triple_fr_emp.csv'

###########################################################
#simulate gene trees under coalescent model for BP species trees

speciesTr_BP_text=readLines(speciesTr_BP_path)
for (i in 1:length(speciesTr_BP_text)){
	geneTr_sim_text=geneTr_sim(speciesTr_BP_text[i],geneTr_num)
}

