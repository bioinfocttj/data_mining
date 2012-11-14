import lxml.etree as etree

file = open('echantillon.xml','r')
tree = etree.parse(file)
root = tree.getroot()

#------recup mots cles--------
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
	
print kw




# -------calul nb alpha beta----------
nb_strand = 0	
nb_turn = 0	
nb_helix = 0
for parent in tree.getiterator():
	print parent
	nb_strand = 0	
	nb_turn = 0	
	nb_helix = 0
	if parent.tag=='entry':
		for child in parent:
			if child.tag=='feature' and child.attrib['type']=='strand':
				nb_strand+=1
			if child.tag=='feature' and child.attrib['type']=='turn':
				nb_turn+=1
			if child.tag=='feature' and child.attrib['type']=='helix':
				nb_helix+=1
		strand = etree.Element("strand")
		strand.attrib["nb"] = str(nb_strand)
		parent.append(strand)
		helix = etree.Element("helix")
		helix.attrib["nb"] = str(nb_helix)
		parent.append(helix)
		turn = etree.Element("turn")
		turn.attrib["nb"] = str(nb_turn)
		parent.append(turn)

# ------supprime les mots cles non selectionnes-----		
for node in tree.getiterator():
	if node.tag=='keyword':
		if node.text not in kw :
			node.getparent().remove(node)
tree.write("Test.xml")
file.close()
