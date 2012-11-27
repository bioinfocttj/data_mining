import lxml.etree as etree #parse XML
import math
import sys
import copy

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

def ajout_identifiant(cluster,matrix):
	ligne = 1
	for key in range(len(cluster)):#ajout des identifiants ds la 1ere ligne et 1ere colonne de la matrice de distance
		matrix[0][ligne]=cluster[key]
		matrix[ligne][0]=cluster[key]
		ligne+=1
	return matrix

def remplir_matrice_initiale_length(dico_prot,matrix):
	for i in range (1,len(matrix)):
		key1=matrix[0][i]
		a = dico_prot[key1[0]]['length']
		for j in range (1,len(matrix)):
			key2=matrix[j][0]
			b = dico_prot[key2[0]]['length']
			matrix[i][j] = dist_length(a,b)
	return matrix

def mini(matrix): #plus petite distance dans la matrice
	savI=1
	savJ=2
	x=2
	mini=sys.maxint
	for i in range (x,len(matrix)):
		for j in range (1,len(matrix)):
			if matrix[i][j]<mini and i!=j:
				mini=matrix[i][j]
				savI=matrix[0][i]
				savJ=matrix[j][0]
		x+=1
	return savI[0],savJ[0]
	
def trouver_distance(matrice,a,b):
	for i in range (1,len(matrice)):
		if a in matrice[0][i]:
			for j in range (1,len(matrice)):
				if b in matrice[j][0]:
					d = matrice[i][j]
	return d

def remplir_matrice(matrice,i,j,cluster): #fusion valeurs i et j 
	newcluster=new_cluster(cluster,i,j)
	matrix = init_matrix(len(matrice)-1)
	matrix = ajout_identifiant(newcluster,matrix)
#	print 'i,j', i, j
	for a in range (1,len(matrix)):
		if i in matrix[0][a]  :
			for b in range(1,len(matrix)):
				if matrix[a][0][0]==matrix[0][b][0]:
					matrix[a][b]=0.0
				else :
					matrix[a][b]=(trouver_distance(matrice,i,matrix[0][b][0])+trouver_distance(matrice,j,matrix[0][b][0]))/2
					matrix[b][a]=matrix[a][b]
		else :
			for b in range(1,len(matrix)):
				if matrix[a][0][0]==matrix[0][b][0]:
					matrix[a][b]=0.0
				else :
					matrix[a][b]=trouver_distance(matrice,matrix[a][0][0],matrix[0][b][0])
					matrix[b][a]=matrix[a][b]
	return matrix,newcluster

def liste_cluster(listeCluster,cluster):
	listeCluster.append(cluster)
	return listeCluster

def creer_cluster(dico_prot):
	cluster=[]
	for key in dico_prot.keys():
		temp=[key]
		cluster.append(temp)
	return cluster

def new_cluster(cluster,i,j):
	newCluster=[]
	for k in cluster :
		if i not in k and j not in k:
			newCluster.append(k)
		if i in k :
			l1 = k
		if j in k :
			l2 = k
	for a in l2:
		l1.append(a)
	newCluster.append(l1)
	return newCluster

def distance_structure(helixA,helixB,strandA,strandB,turnA,turnB):
	d=abs(helixB-helixA)+abs(strandB-strandA)+abs(turnB-turnA)
	return d 

def remplir_matrice_initiale_structure(dico,matrix):
	for i in range (1,len(matrix)):
		key1=matrix[0][i]
		Ahelix = dico[key1[0]]['pct_helix']
		Astrand = dico[key1[0]]['pct_strand']
		Aturn = dico[key1[0]]['pct_turn']
		for j in range (1,len(matrix)):
			key2=matrix[j][0]
			Bhelix = dico[key2[0]]['pct_helix']
			Bstrand = dico[key2[0]]['pct_strand']
			Bturn = dico[key2[0]]['pct_turn']
			d=distance_structure(Ahelix,Bhelix,Astrand,Bstrand,Aturn,Bturn)
			#print 'd=',d
			matrix[i][j] = d
	return matrix

def cluster_taille(dico_prot):
	cluster = []
	for key in dico_prot :
		if dico_prot[key]['length']<100 :
			
	
#-------------------------------------------------------------------------------------------------------------------
#------------------------------------------------                  -------------------------------------------------
#------------------------------------------------       MAIN       -------------------------------------------------
#------------------------------------------------                  -------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
dico_prot,dico_struct,dico_chimie = create_table('mini_echantillon.xml')

#---------- Calcul matrice de distance length + cluster -----------------

print '--------------cluster niveau taille de la chaine----------------'
matrix_length = init_matrix(len(dico_prot)+1)

#---- remplissage de la matrice ----
clust=creer_cluster(dico_prot)
mat=ajout_identifiant(clust,matrix_length)
mat=remplir_matrice_initiale_length(dico_prot,mat)
#for l in mat :
#	print l
print len(clust)#sert a verifier a quel niveaux de cluster on est
clustersTotaux=[]

while (len(clust))>1:
	savI,savJ=mini(mat)
	mat,clust=remplir_matrice(mat,savI,savJ,clust)
	temp=[]
	
	temp = copy.deepcopy(clust)
	#print temp
	print len(clust)
	clustersTotaux.append(temp)
	
print "---------------CLUSTERS TOTAUX LENGTH----------------"
for l in clustersTotaux:
	print l

print '--------------cluster niveau structure secondaires----------------'
matrix_struct = init_matrix(len(dico_struct)+1)
clust2=creer_cluster(dico_struct)
mat2=ajout_identifiant(clust2,matrix_struct)
for l in mat2 :
	print l
mat2=remplir_matrice_initiale_structure(dico_struct,mat2)
print 'matrice remplie'
for l in mat2 :
	print l
print len(clust)#sert a verifier a quel niveaux de cluster on est
clustersTotauxStructure=[]

while (len(clust2))>1:
	savI,savJ=mini(mat2)
	mat2,clust2=remplir_matrice(mat2,savI,savJ,clust2)
	temp=[]
	
	temp = copy.deepcopy(clust2)
	#print temp
	print len(clust2)
	clustersTotauxStructure.append(temp)
print "---------------CLUSTERS TOTAUX STRUCTURE----------------"
for l in clustersTotauxStructure:
	print l

print '--------------cluster niveau structure secondaires----------------'
#print distance_structure(dico_struct['Q15173']['pct_helix'],dico_struct['Q15173']['pct_strand'],dico_struct['Q15173']['pct_turn'],dico_struct['P04229']['pct_helix'],dico_struct['P04229']['pct_strand'],dico_struct['P04229']['pct_turn'])
