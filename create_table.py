import lxml.etree as etree #parse XML
import math
from scipy import cluster
from matplotlib.pyplot import show

file = open('echcalc.xml','r')

tree = etree.parse(file)

dico_prot = {} # dictionnaire correspondant a la table des proteines
dico_struct = {} # dictionnaire correspondant a la table des structures
dico_chimie = {} # dictionnaire correspondant a la table des proprietes chimiques

for parent in tree.findall('.//entry'):
	nb_gene = 0
	list_tissue=[]
	for child in parent :
		if child.tag == 'accession':
			id_prot = child.text
			dico_prot[id_prot]=[0,0,0,0]
			dico_struct[id_prot]=[0,0,0,0,0,0,0]
			dico_chimie[id_prot]=[0,0]
		if child.tag=='fullName':
			dico_prot[id_prot][0] = child.text
		if child.tag=='tissue':
			list_tissue.append(child.text)
		if child.tag =='sequence':
			dico_prot[id_prot][2]=float(child.attrib['length'])
			dico_chimie[id_prot][0]=float(child.attrib['gravy'])
			dico_chimie[id_prot][1]=float(child.attrib['phi'])
			dico_struct[id_prot][0]=float(child.attrib['cysteine'])
		if child.tag == 'gene':
			nb_gene+=1
			if len(child.getchildren())>0:
				dico_prot[id_prot][1]=child.getchildren()[0].text
			else :
				dico_prot[id_prot][1]='NOGENE'	
		if child.tag == 'strand' :
			dico_struct[id_prot][1]=float(child.attrib['nb'])
			dico_struct[id_prot][2]=float(child.attrib['pct'])
		if child.tag == 'helix' :
			dico_struct[id_prot][3]=float(child.attrib['nb'])
			dico_struct[id_prot][4]=float(child.attrib['pct'])
		if child.tag == 'turn' :
			dico_struct[id_prot][5]=float(child.attrib['nb'])
			dico_struct[id_prot][6]=float(child.attrib['pct'])
	if nb_gene == 0 :
		dico_prot[id_prot][1]='NOGENE'
	dico_prot[id_prot][3]=list_tissue



# dico_prot = {id_prot : [fullname,gene,length,[list tissue]]}
# dico_struct = {id_prot : [nb_cystein,nb_strand,pct_strand,nb_helix,pct_helix,nb_turn,pct_turn]}
# dico_chimie = {id_prot : [gravy,phi]}

def init_matrix(d):   #creation matrice de distance de dimension d
	matrix=[]
	for i in range(d):
		line=[]
		for j in range(d):
			line.append(0.0)
		matrix.append(line)
	return matrix
 
def reduct_matrix(matrix):
    # ramene toutes les valeurs d'une matrice entre 0 et 1,
    # permet de faciliter la ponderation des matrices
    red_matrix=init_matrix(len(matrix))
    maxi=0
    for i in range(len(matrix)):
        for j in range(i+1,len(matrix)):
            if matrix[i][j]>maxi:
                maxi=matrix[i][j]
    if maxi!=0:
        for i in range(len(matrix)):
            for j in range(i,len(matrix)):
                red_dist=(matrix[i][j])/maxi
                red_matrix[i][j]=red_dist
                red_matrix[j][i]=red_dist
    else:
        red_matrix=matrix
    return red_matrix 
    
def dist_length(a,b): # calcul distance 
	dist = math.fabs(a-b)
	return dist

matrix_length = init_matrix(len(dico_prot))
list_id = []  #pour garder en memoire l'ordre des proteines dans la matrice

for key in dico_prot :
	list_id.append(key)
i=0
j=0		
for i in range (len(matrix_length)):
	a = dico_prot[list_id[i]][2]
	for j in range (len(matrix_length)):
		b = dico_prot[list_id[j]][2]
		matrix_length[i][j] = dist_length(a,b)
		
matrix_red = reduct_matrix(matrix_length)	

print matrix_red

