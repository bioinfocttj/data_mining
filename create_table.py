import lxml.etree as etree #parse XML
import math
#from scipy import cluster
#from matplotlib.pyplot import show

def create_table(file):

	file = open(file,'r')

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
	return dico_prot,dico_struct,dico_chimie


# dico_prot = {id_prot : [fullname,gene,length,[list tissue]]}
# dico_struct = {id_prot : [nb_cystein,nb_strand,pct_strand,nb_helix,pct_helix,nb_turn,pct_turn]}
# dico_chimie = {id_prot : [gravy,phi]}



dico_prot,dico_struct,dico_chimie = create_table('echantillon_final.xml')

#---------- Calcul matrice de distance length + cluster -----------------

def init_matrix(d):   #creation matrice de distance de dimension d
	matrix=[]
	for i in range(d):
		line=[]
		for j in range(d):
			line.append(0.0)
		matrix.append(line)
	return matrix

def dist_length(a,b): # calcul distance 
	dist = math.fabs(a-b)
	return dist


matrix_length = init_matrix(len(dico_prot)+1)

#---- remplissage de la matrice ----

ligne = 1
for key in dico_prot.keys():#ajout des identifiants ds la 1ere ligne et 1ere colonne de la matrice de distance
	matrix_length[0][ligne]=key
	matrix_length[ligne][0]=key
	ligne+=1

for i in range (1,len(matrix_length)):
	key1=matrix_length[0][i]
	a = dico_prot[key1]['length']
	for j in range (1,len(matrix_length)):
		key2=matrix_length[j][0]
		b = dico_prot[key2]['length']
		matrix_length[i][j] = dist_length(a,b)

def min(matrix_length):
	savI=1
	savJ=2
	x=2
	mini=matrix_length[1][2]
	for i in range (x,len(matrix_length)):
		for j in range (1,len(matrix_length)):
			if matrix_length[i][j]<mini:
				mini=matrix_length[i][j]
				savI=i
				savJ=j
		x+=1
	return savI,savJ

def matrice_distance(matrice,i,j):
	matrix = init_matrix(len(matrice)-1)
	for ligne in range (len(matrice)):
		matrix[][]=(matrice[i][ligne]+matrice[j][ligne])/2
