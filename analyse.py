import lxml.etree as etree #parse XML
import numpy
from mpl_toolkits.mplot3d import Axes3D

from pylab import *
 

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
	
def stat(file_name,dico,texte):
	file = open(file_name,'r')
	for line in file :
		if line!=('\n'):
			valeurs=[]
			tab = line.split(' ')
			for i in range (len(tab)-1):
				valeurs.append(dico[tab[i]][texte])
			moyenne = numpy.average(valeurs)
			ecartType = numpy.std(valeurs)
			print 'moyenne', numpy.average(valeurs)
			print 'ecart-type',numpy.std(valeurs)
	file.close()

dico_prot,dico_struct,dico_chimie = create_table('homo_final.xml')

def resultat(file_name,file_result):
	nb_cluster_seul = 0
	nb_cluster=0
	file = 	file = open(file_name,'r')
	file2 = open(file_result,'w')
	for line in file :
		if line!=('\n'):
			nb_cluster+=1
			s = 'cluster n'+str(nb_cluster)+('\n')
			file2.write(s)
			tab = line.split(' ')
			if (len(tab))==2:
				nb_cluster_seul+=1
			for i in range (len(tab)-1):
				s = '	-'+str(dico_prot[tab[i]]['fullname'])+('\n')
				file2.write(s)
	file2.write('Nombre de cluster a une proteine')
	file2.write(str(nb_cluster_seul))
	print nb_cluster
	file.close()
	file2.close()

def moyenne_intracluster(file_name,file_result):
	nb_cluster=0
	file = 	file = open(file_name,'r')
	file2 = open(file_result,'w')
	hydro_std= []
	histo_hydro = []
	for line in file :
		taille = []
		hydro = []
		phi = []
		cys = []
		turn = []
		helix = []
		strand = []
		if line!=('\n'):
			nb_cluster+=1
			s = 'cluster n '+str(nb_cluster)+('\n')
			file2.write(s)
			tab = line.split(' ')
			for i in range (len(tab)-1):
				helix.append(dico_struct[tab[i]]['pct_helix'])
				turn.append(dico_struct[tab[i]]['pct_turn'])
				strand.append(dico_struct[tab[i]]['pct_strand'])
				taille.append(dico_prot[tab[i]]['length'])
				hydro.append(dico_chimie[tab[i]]['gravy'])
				phi.append(dico_chimie[tab[i]]['phi'])
				cys.append(dico_struct[tab[i]]['nb_cystein'])
			moyenne_struct = numpy.average(taille)
			moyenne_taille = numpy.average(taille)
			moyenne_hydro = numpy.average(hydro)
			moyenne_phi = numpy.average(phi)
			moyenne_cys = numpy.average(cys) 
			moyenne_strand = numpy.average(strand)
			moyenne_turn = numpy.average(turn)
			moyenne_helix = numpy.average(helix)
			std_struct = numpy.std(taille)
			std_taille = numpy.std(taille)
			std_hydro = numpy.std(hydro)
			std_phi = numpy.std(phi)
			std_cys = numpy.std(cys) 
			std_strand = numpy.std(strand)
			std_turn = numpy.std(turn)
			std_helix = numpy.std(helix)
			file2.write('moyenne taille :')
			file2.write(str(moyenne_taille))
			file2.write('ecart-type taille :')
			file2.write(str(std_taille))
			file2.write('\n')
			file2.write('moyenne hydro :')
			file2.write(str(moyenne_hydro))
			file2.write('ecart-type taille :')
			file2.write(str(std_hydro))
			file2.write('\n')
			file2.write('moyenne phi :')
			file2.write(str(moyenne_phi))
			file2.write(' ecart-type taille :')
			file2.write(str(std_phi))
			file2.write('\n')
			file2.write('moyenne cysteine : ')
			file2.write(str(moyenne_cys))
			file2.write(' ecart-type taille : ')
			file2.write(str(std_cys))
			file2.write('\n')
			file2.write('moyenne beta : ')
			file2.write(str(moyenne_helix))
			file2.write(' ecart-type taille : ')
			file2.write(str(std_helix))
			file2.write('\n')
			file2.write('moyenne alpha :')
			file2.write(str(moyenne_strand))
			file2.write(' ecart-type taille : ')
			file2.write(str(std_strand))
			file2.write('\n')
			file2.write('moyenne coude :')
			file2.write(str(moyenne_turn))
			file2.write(' ecart-type taille : ')
			file2.write(str(std_turn))
			file2.write('\n')
			hydro_std.append(std_hydro)
			histo_hydro.append(moyenne_hydro)
	file.close()
	file2.close()
	y=[]
	for i in range(len(histo_taille)):
		y.append(i)
	fig = figure()
	ax = fig.add_subplot(1,1,1)
	ax.bar(y,histo_hydro,yerr=hydro_std)
	xlabel('clusters')
	ylabel('hydrophobicite')
	title("moyenne et ecart-type intercluster pour l'hydrophobicite")
	#~ show()#---------------- decommenter pour afficher les histogrammes

moyenne_intracluster('resultat_taille1.txt','resultat_intracluster1.txt')
moyenne_intracluster('resultat_structure2.txt','resultat_intracluster2.txt')
moyenne_intracluster('resultat_hydro3.txt','resultat_intracluster3.txt')
moyenne_intracluster('resultat_phi4.txt','resultat_intracluster4.txt')
moyenne_intracluster('resultat_cys5.txt','resultat_intracluster5.txt')
moyenne_intracluster('resultat_organe6.txt','resultat_intracluster6.txt')


