import lxml.etree as etree #parse XML

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


	



