
file = open('/net/cremi/alelievr/espaces/travail/data_mining/echantillon.xml','r')
file2 = open('echantillon_parse.xml','w')

file2.write('<uniprot>')
for line in file :
	if '<entry' in line :
		file2.write(line) #pour separer les proteines
	if '</entry>' in line :
		file2.write(line) 
	if '<accession>' in line :
		file2.write(line)
	if '<name>' in line :
		file2.write(line)
	if '<fullName' in line :
		file2.write(line)
	if '<tissue' in line :
		file2.write(line)
	if '<keyword' in line :
		file2.write(line)
	if '<sequence length' in line :
		while '</sequence>' not in line:
			file2.write(line)
			line = file.next()
		file2.write(line)
	if '<subcellularLocation>' in line:
		line = file.next()
		if '<location>' in line:
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
