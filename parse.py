
file = open('/net/cremi/alelievr/test.xml','r')
file2 = open('test','w')

for line in file :
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
	if '<sequence' in line :
		line = file.next()
		while line!='</sequence>':
			file2.write(line)
			line = file.next()

file.close()
file2.close()
