import lxml.etree as etree

file = open('echantillon_parse.xml','r')
tree = etree.parse(file)
root = tree.getroot()


# -------calul nb alpha beta----------
nb_strand = 0	
nb_turn = 0	
nb_helix = 0
for parent in tree.getiterator():
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

# ------suppression tissue appaissant plusieurs fois pour une meme proteine-----	
	
for parent in tree.getiterator():
    if parent.tag=='entry':
        tissue=[]
        for child in parent:
            if child.tag=='tissue':
                if child.text not in tissue :
                    tissue.append(child.text)
                else:
                    child.getparent().remove(child)
			

tree.write("echcalc.xml")

