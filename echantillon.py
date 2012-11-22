
file = open('/net/stockage/bioinfo_promo_2011_2013/DataMining/homo_sapiens.xml','r')
file2 = open('echantillon.xml','w')

cpt=0
cpt2=0
file2.write('<uniprot>')
for line in file :
	if '<entry' in line :
		if cpt%10==0:
			while '</entry>' not in line:
				file2.write(line)
				line = file.next()
			file2.write(line)
			cpt2+=1
		cpt+=1
		
file2.write('</uniprot>')		
print "CPT :" , cpt2
print "CPT prot :" , cpt
