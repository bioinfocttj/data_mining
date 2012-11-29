import lxml.etree as etree #parse XML
import math
import sys
import copy
import os 

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






#---------- Calcul matrice de distance length + cluster -----------------


def init_matrix(d):   #creation matrice de distance de dimension d
	matrix=[]
	for i in range(d):
		line=[]
		for j in range(d):
			line.append(0.0)
		matrix.append(line)
	return matrix

def dist(a,b): # calcul distance 
	dist = math.fabs(a-b)
	return dist

#~ matrix_length = init_matrix(len(dico_prot)+1)

#---- remplissage de la matrice ----

def tri_taille(dico): #premier etape de clustarisation (taille) taille de cluster redefinies arbitrairement
	liste=[[],[]]
	for key in dico :
		liste[0].append(key)
		liste[1].append(dico[key]['length'])
	liste=Tri_Bulle(liste)
	return liste
	
def Tri_Bulle(liste):
	l=liste[1]
	l1=liste[0]
	for i in xrange(len(l)):
		for j in reversed(xrange(i,len(l))):
			if l[j]<l[j-1]:
				t=l[j]
				t1=l1[j]
				l[j]=l[j-1]
				l1[j]=l1[j-1]
				l[j-1]=t
				l1[j-1]=t1
	return liste
	


def ajout_identifiant(cluster,matrix):
	ligne = 1
	for key in range(len(cluster)):#ajout des identifiants ds la 1ere ligne et 1ere colonne de la matrice de distance
		matrix[0][ligne]=cluster[key]
		matrix[ligne][0]=cluster[key]
		ligne+=1
	return matrix

def remplir_matrice_cluster(dico,matrix,texte):
	for i in range (1,len(matrix)):
		key1=matrix[0][i]
		a = dico[key1[0]][texte]
		for j in range (1,len(matrix)):
			key2=matrix[j][0]
			b = dico[key2[0]][texte]
			matrix[i][j] = dist(a,b)
	return matrix

def remplir_matrice_organe(dico,matrix,texte):
	for i in range (1,len(matrix)):
		key1=matrix[0][i]
		a = dico[key1[0]][texte]
		for j in range (1,len(matrix)):
			key2=matrix[j][0]
			b = dico[key2[0]][texte]
			matrix[i][j] = distance_liste(a,b)
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

def creer_cluster(liste):
	cluster=[]
	for i in liste:
		temp=[i]
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

def distance_liste(l1,l2):
	intersection = list(set(l1) & set(l2))
	if (len(l1))>=(len(l2)):
		L = len(l1)
	else : 
		L = len(l2)
	if L!=0 :
		return (len(intersection)/L)
	else :
		return 0


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
			matrix[i][j] = d
	return matrix

def clustering_taille(liste):
	cluster=[]
	temp=[]
	for i in range (len(liste[0])):
		temp.append(liste[0][i])
		if i%200==0 and i!=0:
			cluster.append(temp)
			temp=[]
	return cluster
	
def cluster_struct(liste):
	matrix = init_matrix(len(liste)+1)
	cluster=creer_cluster(liste)
	matrix = ajout_identifiant(cluster,matrix)
	matrix=remplir_matrice_initiale_structure(dico_struct,matrix)
	temp = []
	while (len(cluster))>4:
		print len(cluster)
		savI,savJ=mini(matrix)
		matrix,cluster=remplir_matrice(matrix,savI,savJ,cluster)
		temp=[]
		temp = copy.deepcopy(cluster)
	return temp
	
def cluster_autre(liste,dico,texte):
	matrix = init_matrix(len(liste)+1)
	cluster=creer_cluster(liste)
	matrix = ajout_identifiant(cluster,matrix)
	matrix=remplir_matrice_cluster(dico,matrix,texte)
	temp = []
	if (len(cluster))<4 :
		temp = copy.deepcopy(cluster)
	while (len(cluster))>=4:
		print len(cluster)
		savI,savJ=mini(matrix)
		matrix,cluster=remplir_matrice(matrix,savI,savJ,cluster)
		temp=[]
		temp = copy.deepcopy(cluster)
	return temp

def cluster_organe(liste,dico,texte):
	matrix = init_matrix(len(liste)+1)
	cluster=creer_cluster(liste)
	matrix = ajout_identifiant(cluster,matrix)
	matrix=remplir_matrice_organe(dico,matrix,texte)
	temp = []
	if (len(cluster))<4 :
		temp = copy.deepcopy(cluster)
	while (len(cluster))>4:
		print len(cluster)
		savI,savJ=mini(matrix)
		matrix,cluster=remplir_matrice(matrix,savI,savJ,cluster)
		temp=[]
		temp = copy.deepcopy(cluster)
	return temp
	
def ecriture(cluster,fichier):
	file = open(fichier,'w')
	if cluster==clusterTaille:
		for i in range(len(cluster)):
			file.write('\n')
			for j in range (len(cluster[i])):
				file.write(cluster[i][j])
				file.write(' ')
	elif cluster==clusterStruct:
		for i in range(len(cluster)):
			file.write('\n')
			for j in range (len(cluster[i])):
				file.write('\n')
				for k in range (len(cluster[i][j])):
					file.write(cluster[i][j][k])
					file.write(' ')
	elif cluster==clusterHydro:
		for i in range(len(cluster)):
			file.write('\n')
			for j in range (len(cluster[i])):
				file.write('\n')
				for k in range (len(cluster[i][j])):
					file.write('\n')
					for l in range (len(cluster[i][j][k])):
						file.write(cluster[i][j][k][l])
						file.write(' ')
	elif cluster==clusterPhi:
		for i in range(len(cluster)):
			file.write('\n')
			for j in range (len(cluster[i])):
				file.write('\n')
				for k in range (len(cluster[i][j])):
					file.write('\n')
					for l in range (len(cluster[i][j][k])):
						file.write('\n')
						for m in range (len(cluster[i][j][k][l])):
							file.write(cluster[i][j][k][l][m])
							file.write(' ')
	elif cluster==clusterCys:
		for i in range(len(cluster)):
			file.write('\n')
			for j in range (len(cluster[i])):
				file.write('\n')
				for k in range (len(cluster[i][j])):
					file.write('\n')
					for l in range (len(cluster[i][j][k])):
						file.write('\n')
						for m in range (len(cluster[i][j][k][l])):
							file.write('\n')
							for n in range (len(cluster[i][j][k][l][m])):
								file.write(cluster[i][j][k][l][m][n])
								file.write(' ')
	elif cluster==clusterOrgane:
		for i in range(len(cluster)):
			file.write('\n')
			for j in range (len(cluster[i])):
				file.write('\n')
				for k in range (len(cluster[i][j])):
					file.write('\n')
					for l in range (len(cluster[i][j][k])):
						file.write('\n')
						for m in range (len(cluster[i][j][k][l])):
							file.write('\n')
							for n in range (len(cluster[i][j][k][l][m])):
								file.write('\n')
								for o in range (len(cluster[i][j][k][l][m][n])):
									file.write(cluster[i][j][k][l][m][n][o])
									file.write(' ')
		
	file.close()
	

#-------------------------------------------------------------------------------------------------------------------
#------------------------------------------------                  -------------------------------------------------
#------------------------------------------------       MAIN       -------------------------------------------------
#------------------------------------------------                  -------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

dico_prot,dico_struct,dico_chimie = create_table('echantillon_final.xml')

# -------- Clusterisation Taille -------------------
l = tri_taille(dico_prot)
clusterTaille = clustering_taille(l)    #----- clusterTaille = clusters de taille
clusterStruct = [] #----- clusterStruct = clusters de structure a partir des clusters de taille
 #----- clusterHydro = clusters d'hydrophobicite a partir des clusters de structure
clusterPhi = [] #----- clusterPhi = clusters Phi a partir des clusters de d'hydrophobicite
clusterCys = []
clusterOrgane = []
# -------- Clusterisation Structure-------------------
for i in range (len(clusterTaille)) :
	clusterStruct.append(cluster_struct(clusterTaille[i]))


#~ # -------- Clusterisation Hydrophobicite------------------
clusterHydro = copy.deepcopy(clusterStruct)
for i in range (len(clusterStruct)) :
	for j in range (len(clusterStruct[i])):
		clusterHydro[i][j]=(cluster_autre(clusterStruct[i][j],dico_chimie,'gravy'))



# -------------- Clusterisation Phi ----------------------
clusterPhi = copy.deepcopy(clusterHydro)
for i in range (len(clusterHydro)) :
	for j in range (len(clusterHydro[i])):
		for k in range (len(clusterHydro[i][j])):
			clusterPhi[i][j][k]=(cluster_autre(clusterHydro[i][j][k],dico_chimie,'phi'))
		
		
#-----------------Clusterisation Cysteine---------------------
clusterCys = copy.deepcopy(clusterPhi)
for i in range (len(clusterPhi)) :
	for j in range (len(clusterPhi[i])):
		for k in range (len(clusterPhi[i][j])):
			for l in range (len(clusterPhi[i][j][k])):
				clusterCys[i][j][k][l]=(cluster_autre(clusterPhi[i][j][k][l],dico_struct,'nb_cystein'))
				
#~ #-----------------Clusterisation organe---------------------
clusterOrgane = copy.deepcopy(clusterCys)
for i in range (len(clusterCys)) :
	for j in range (len(clusterCys[i])):
		for k in range (len(clusterCys[i][j])):
			for l in range (len(clusterCys[i][j][k])):
				for m in range (len(clusterCys[i][j][k][l])):
					clusterOrgane[i][j][k][l][m]=(cluster_organe(clusterCys[i][j][k][l][m],dico_prot,'tissue')) 

#----------------Ecriture resultats fichiers -----------------

ecriture(clusterTaille,'resultat_taille1.txt')
ecriture(clusterStruct,'resultat_structure2.txt')
ecriture(clusterHydro,'resultat_hydro3.txt')
ecriture(clusterPhi,'resultat_phi4.txt')
ecriture(clusterCys,'resultat_cys5.txt')
ecriture(clusterOrgane,'resultat_organe6.txt')

#~ print '-----------------------------------------------------------clusterTaille'
#~ print clusterTaille
#~ 
#~ print '-----------------------------------------------------------clusterStruct'
#~ print clusterStruct
#~ 
#~ 
#~ print '-----------------------------------------------------------clusterHydro'
#~ print clusterHydro
#~ 
#~ print '-----------------------------------------------------------clusterPhi'
#~ print clusterPhi
#~ 
#~ print '-----------------------------------------------------------clusterCys'
#~ print clusterCys
#~ 
#~ 
#~ print '-----------------------------------------------------------clusterOrgane'
#~ print clusterOrgane



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
for i in range (len(clusterTaille)):
	file.write('Proteines')
	file.write('-> ')
	n = conc_liste(clusterTaille[i])
	file.write(n)
	file.write('\n')
	file.write(n)
	if len(clusterTaille[i])!=1:
		file.write('[shape=point]')
	else :
		file.write('[shape=point]')		
	file.write('\n')
	for j in range (len(clusterStruct[i])) :
		n1 = conc_liste(clusterStruct[i][j])
		file.write(n)
		file.write('-> ')
		file.write(n1)
		file.write('\n')
		file.write(n1)
		if len(clusterStruct[i][j])!=1:
			file.write('[shape=point]')
		else :
			file.write('[shape=point]')		
		file.write('\n')
		for k in range (len(clusterHydro[i][j])) :
			n2 = conc_liste(clusterHydro[i][j][k])
			if len(clusterHydro[i][j][k])!=1:
				file.write('[shape=point]')
			if n1!=n2 :
				file.write(n1)
				file.write('-> ')
				file.write(n2)
				file.write('\n')
				file.write(n2)
			else :
				file.write('[shape=point]')	
			file.write('\n')
			for l in range (len(clusterPhi[i][j][k])) :
				n3 = conc_liste(clusterPhi[i][j][k][l])
				if len(clusterPhi[i][j][k][l])!=1:
					file.write('[shape=point]')
				if n2!=n3 :
					file.write(n2)
					file.write('-> ')
					file.write(n3)
					file.write('\n')
					file.write(n3)
				
				else :
					file.write('[shape=point]')	
				file.write('\n')

file.write('}')


file.close()

os.popen('dotty resultat.dot&')
