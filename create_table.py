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
			dico_prot[id_prot]={'fullname':0,'gene':0,'length':0,'tissue':0}
			dico_struct[id_prot]={'nb_cystein':0,'nb_strand':0,'pct_strand':0,'nb_helix':0,'pct_helix':0,'nb_turn':0,'pct_turn':0}
			dico_chimie[id_prot]={'gravy':0,'phi':0}
		if child.tag=='fullName':
			dico_prot[id_prot]['fullname'] = child.text
		if child.tag=='tissue':
			list_tissue.append(child.text)
		if child.tag =='sequence':
			dico_prot[id_prot]['length']=float(child.attrib['length'])
			dico_chimie[id_prot]['gravy']=float(child.attrib['gravy'])
			dico_chimie[id_prot]['phi']=float(child.attrib['phi'])
			dico_struct[id_prot]['nb_cystein']=float(child.attrib['cysteine'])
		if child.tag == 'gene':
			nb_gene+=1
			if len(child.getchildren())>0:
				dico_prot[id_prot]['gene']=child.getchildren()[0].text
			else :
				dico_prot[id_prot]['gene']='NOGENE'
		if child.tag == 'strand' :
			dico_struct[id_prot]['nb_strand']=float(child.attrib['nb'])
			dico_struct[id_prot]['pct_strand']=float(child.attrib['pct'])
		if child.tag == 'helix' :
			dico_struct[id_prot]['nb_helix']=float(child.attrib['nb'])
			dico_struct[id_prot]['pct_helix']=float(child.attrib['pct'])
		if child.tag == 'turn' :
			dico_struct[id_prot]['nb_turn']=float(child.attrib['nb'])
			dico_struct[id_prot]['pct_turn']=float(child.attrib['pct'])
	if nb_gene == 0 :
		dico_prot[id_prot]['gene']='NOGENE'
	dico_prot[id_prot]['tissue']=list_tissue



# dico_prot = {id_prot : [fullname,gene,length,[list tissue]]}
# dico_struct = {id_prot : [nb_cystein,nb_strand,pct_strand,nb_helix,pct_helix,nb_turn,pct_turn]}
# dico_chimie = {id_prot : [gravy,phi]}





#---------- Calcul matrice de distance length + cluster -----------------

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

#---- remplissage de la matrice ----
i=0
j=0
for key1 in dico_prot.keys():
	a = dico_prot[key1]['length']
	for key2 in dico_prot.keys():
		b = dico_prot[key2]['length']
		matrix_length[i][j] = dist_length(a,b)
		j+=1
	i+=1
	j=0
	


#matrix_red = reduct_matrix(matrix_length)	

#---- tests ----
"""
step_cluster =[] #contient toutes les etapes de clusterisation
cluster = [] # cluster en cours

# toutes les proteines sont dans un seul cluster
for key in dic_prot :
	cluster.append(key)
	
step_cluster.append(cluster) #ajout premiere etape
i=0
j=0

while  (len(cluster))!=0:
	d=0
	for key1 in dico_prot.keys():
		for key2 in range dico_prot.keys():
			j+=1
			if matrix_red[i][j]>d
			d = matrix_red[i][j]
			p1_id = dico_prot[key1]
			p2_id = dico_prot[key2]
	i+=1
	j=0
"""


Y =  cluster.hierarchy.linkage(matrix_length,method='complete')

#Coupage de l'arbre -> recuperation des clusters



Z = cluster.hierarchy.dendrogram(Y,p=10,truncate_mode='lastp')
show()
