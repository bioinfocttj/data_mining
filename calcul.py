import lxml.etree as etree #parse XML
from Bio.SeqUtils import ProtParam #calcul phi

file = open('echantillon_parse.xml','r')
tree = etree.parse(file)
root = tree.getroot()

#----suppression numeros accession et full name si plusieurs-----
for parent in tree.getiterator():
	nb_accession=0
	if parent.tag=='entry':
		for child in parent :
			if child.tag=='accession':
				nb_accession+=1
				if nb_accession>1:
					child.getparent().remove(child)
					
for parent in tree.getiterator():
	nb_name=0
	if parent.tag=='entry':
		for child in parent :
			if child.tag=='fullName':
				nb_name+=1
				if nb_name>1:
					child.getparent().remove(child)					
		
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
			
# ------calcul nb cysteine-----

for node in tree.getiterator():
	cpt_cys=0
	if node.tag =='sequence':
		sequence = node.text
		for aa in sequence :
			if aa=='C':
				cpt_cys+=1
		node.attrib["cysteine"] = str(cpt_cys)
	
# -------calcul pourcentage aa hydrophobes ------		

aa_hphobes =  ['G', 'A', 'V','L', 'I', 'F', 'W', 'Y']   # pas sure des mes sources 

for node in tree.getiterator():
	cpt_hphobe=0
	if node.tag =='sequence':
		sequence = node.text
		for aa in sequence :
			if aa in aa_hphobes :
				cpt_hphobe+=1
		length = float(node.attrib.get('length'))
		pourcent_hphobe = cpt_hphobe / length * 100
		node.attrib["hydrophobicity"] = str(round(pourcent_hphobe,1)) #  1 chiffre apres la virgule 
	
# -------------calcul PHi-----------


for node in tree.findall('.//sequence') :
	sequence = node.text  
	params  = ProtParam.ProteinAnalysis(sequence)
	phi = params.isoelectric_point()
	node.attrib["phi"] = str(round(phi,2))
	
	
	
	
	
tree.write("echcalc.xml")

