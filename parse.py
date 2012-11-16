
file = open('/net/cremi/tgauci/Documents/datamining/data_mining/echantillon.xml','r')
file2 = open('echantillon_parse.xml','w')

file2.write('<uniprot>')
print 'coucou'
for line in file :
	if '<entry' in line :
		file2.write(line) #pour separer les proteines
	if '</entry>' in line :
		file2.write(line) 
	if '<accession>' in line : #numero d'accession de la proteine
		file2.write(line)
	if '<name>' in line : #nom de la proteine
		file2.write(line)
	if '<fullName' in line : #nom complet de la proteine
		file2.write(line)
	if '<tissue' in line : #loclisation au niveau organe
		file2.write(line)
	if '<sequence length' in line : #donnees sur la sequence d'AA (longueur, enchainement AA...)
		while '</sequence>' not in line:
			file2.write(line)
			line = file.next()
		file2.write(line)
	if '<gene>' in line: #gene primaire codant pour la proteine
		file2.write(line)
		while '</gene>' not in line:
			if 'primary' in line:
				file2.write(line)
			line = file.next()
		file2.write(line)
				
	if '<feature type="strand">' in line :   #structures secondaires
		while '</feature>' not in line:
			file2.write(line)
			line = file.next()
		file2.write(line)
	if '<feature type="helix">' in line :   
		while '</feature>' not in line:
			file2.write(line)
			line = file.next()
		file2.write(line)
	if '<feature type="turn">' in line :   
		while '</feature>' not in line:
			file2.write(line)
			line = file.next()
		file2.write(line)

file2.write('</uniprot>')
file.close()
file2.close()
