
file = open('/home/typhaine/Documents/data_mining/P35443.xml','r')
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
		#line = file.next()
		while '</sequence>' not in line:
			file2.write(line)
			line = file.next()
		file2.write(line)
	if '<subcellularLocation>' in line:
		line = file.next()
		if '<location>' in line:
			print line
			file2.write(line)
file.close()
file2.close()
