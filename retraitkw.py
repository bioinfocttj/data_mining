file = open('test','r') #généré à l'exécution de parse.py
file2 = open('echantillonkw.xml','w')
file3 = open('listekw.txt','r')

kw=[]
for line in file3:
	kw.append(line)
	
for i in range (len(kw)) :
	tab = kw[i].split('\n')
	kw[i]=tab[0]
	
for line in file:
	if '<keyword' not in line:
		file2.write(line)
	if '<keyword' in line:
		for i in range(len(kw)):
			if kw[i] in line:
				file2.write(line)
				
file.close()
file2.close()
file3.close()
