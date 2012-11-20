import lxml.etree as etree #parse XML

file = open('echcalc.xml','r')

tree = etree.parse(file)

dico_prot = {} # dictionnaire correspondant a la table des proteines
dico_struct = {} # dictionnaire correspondant a la table des structures
dico_chimie = {} # dictionnaire correspondant a la table des proprietes chimiques

for parent in tree.findall('.//entry'):
	list_tissue=[]
	for child in parent :
		if child.tag == 'accession':
			id_prot = child.text
			dico_prot[id_prot]=[]
			dico_struct[id_prot]=[]
			dico_chimie[id_prot]=[]
		if child.tag=='fullName':
			dico_prot[id_prot].append(child.text)
		if child.tag=='tissue':
			list_tissue.append(child.text)
		if child.tag =='sequence':
			dico_prot[id_prot].append(float(child.attrib['length']))
			dico_chimie[id_prot].append(float(child.attrib['gravy']))
			dico_chimie[id_prot].append(float(child.attrib['phi']))
			dico_struct[id_prot].append(float(child.attrib['cysteine']))
		if child.tag == 'gene':
			if child.getchildren()!=None:
				for e in child.getchildren():
					dico_prot[id_prot].append(e.text)
			else :
				dico_prot[id_prot].append(None)	
		if child.tag == 'strand' or child.tag == 'helix' or child.tag == 'turn':
			dico_struct[id_prot].append(float(child.attrib['nb']))
			dico_struct[id_prot].append(float(child.attrib['pct']))
	dico_prot[id_prot].append(list_tissue)



# dico_prot = {id_prot : [fullname,gene,length,[list tissue]]}
# dico_struct = {id_prot : [nb_cystein,nb_strand,pct_strand,nb_helix,pct_helix,nb_turn,pct_turn]}
# dico_chimie = {id_prot : [gravy,phi]}
