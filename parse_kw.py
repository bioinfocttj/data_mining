
file = open('echantillon.xml','r')

dico= {}
for line in file :
	if '<keyword' in line :
		if line in dico.keys():
			dico[line]+=1
		else :
			dico[line]=1
			
kw=[]
cpt=0
for key in dico.keys():
	if dico[key]>=40:
		cpt+=1		
		kw.append(key)
print cpt

for i in range (len(kw)) :
	tab = kw[i].split('">')
	kw[i]=tab[1]
	
for i in range (len(kw)) :
	tab = kw[i].split('</')
	kw[i]=tab[0]

	
file.close()
