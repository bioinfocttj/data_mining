import os

Liste = ['P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13','P14','P15','P16','P17','P18','P19','P20']
Taille = [['P1','P2','P3','P4','P5'],['P6','P7','P8','P9','P10'],['P11','P12','P13','P14','P15'],['P16','P17','P18','P19','P20']]
Struct = [[['P1','P2','P3'],['P4','P5']],[['P6','P7','P8'],['P9','P10']],[['P11','P12'],['P13','P14','P15']],[['P16','P17'],['P18','P19','P20']]]
PHI = [[[['P1'],['P2'],['P3']],[['P4'],['P5']]],[[['P6'],['P7','P8']],[['P9'],['P10']]],[[['P11'],['P12']],[['P13','P14'],['P15']]],[[['P16'],['P17']],[['P18'],['P19','P20']]]]
Hydro = [[[[['P1']],[['P2']],[['P3']]],[[['P4']],[['P5']]]],[[[['P6']],[['P7'],['P8']]],[[['P9']],[['P10']]]],[[[['P11']],[['P12']]],[[['P13'],['P14']],[['P15']]]],[[[['P16']],[['P17']]],[[['P18']],[['P19'],['P20']]]]]


def conc_liste(liste):
	s=''
	if (len(liste))==0:
		s=liste[0]
	else :
		for i in range (len(liste)) :
			s+=str(liste[i])
		
	return s

file = open('resultat.dot','w')
file.write('digraph cluster {')
file.write('\n')
file.write('\n')
for i in range (len(Taille)):
	file.write('Proteines')
	file.write('-> ')
	n = conc_liste(Taille[i])
	file.write(n)
	file.write('\n')
	file.write(n)
	if len(Taille[i])!=1:
		file.write('[shape=point]')
	else :
		file.write('[shape=plaintext]')		
	file.write('\n')
	for j in range (len(Struct[i])) :
		n1 = conc_liste(Struct[i][j])
		file.write(n)
		file.write('-> ')
		file.write(n1)
		file.write('\n')
		file.write(n1)
		if len(Struct[i][j])!=1:
			file.write('[shape=point]')
		else :
			file.write('[shape=plaintext]')		
		file.write('\n')
		for k in range (len(PHI[i][j])) :
			n2 = conc_liste(PHI[i][j][k])
			if n1!=n2 :
				file.write(n1)
				file.write('-> ')
				file.write(n2)
				file.write('\n')
				file.write(n2)
			if len(PHI[i][j][k])!=1:
				file.write('[shape=point]')
			else :
				file.write('[shape=plaintext]')	
			file.write('\n')
			for l in range (len(Hydro[i][j][k])) :
				n3 = conc_liste(Hydro[i][j][k][l])
				if n2!=n3 :
					file.write(n2)
					file.write('-> ')
					file.write(n3)
					file.write('\n')
					file.write(n3)
				if len(Hydro[i][j][k][l])!=1:
					file.write('[shape=point]')
				else :
					file.write('[shape=plaintext]')	
				file.write('\n')

file.write('}')


file.close()

os.popen('xdot resultat.dot&')
