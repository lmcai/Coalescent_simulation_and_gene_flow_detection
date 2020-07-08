import sys,os
import fnmatch

print('Comparing triplet frequency in the empirical gene trees and simulated gene trees\n')
#store triplet frequency of the simulated dataset in dict
filename=[]
for fn in os.listdir('geneTr_sim'):
    if fnmatch.fnmatch(fn,'*.tsv'):filename.append(fn)

BPtrp_dict={}
for file in filename:
	BPtrp=open('geneTr_sim/'+file).readlines()
	for l in BPtrp:
		a=[int(j) for j in l.split()[1:]]
		a.sort(reverse=True)
		try:
			BPtrp_dict[l.split()[0]].append(a[-2]-a[-1])
		except KeyError:
			BPtrp_dict[l.split()[0]]=[a[-2]-a[-1]]

print('Writing unbalanced triplets to unbalanced.trp.tsv\n')
out=open('unbalanced.trp.tsv','a')

MLtrp=open(sys.argv[1].split('.')[0]+'.trp.tsv').readlines()
MLtrp_dict={}
#unbalanced=[]
for l in MLtrp:
	a=[int(j) for j in l.split()[1:]]
	a.sort(reverse=True)
	if a[-2]-a[-1]>max(BPtrp_dict[l.split()[0]]):
		out.write('\t'.join(l.split())+'\t'+str(max(BPtrp_dict[l.split()[0]]))+'\n')
		#unbalanced.append(l.split()[0])

out.close()
	



