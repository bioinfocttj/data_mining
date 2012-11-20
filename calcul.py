import lxml.etree as etree #parse XML
import sys
sys.path.append('/net/cremi/tgauci/Documents/datamining/data_mining/biopython160/')
from Bio.SeqUtils import ProtParam

file = open('homo_parse.xml','r')
tree = etree.parse(file)
root = tree.getroot()

# -------calcul pourcentage aa hydrophobes ------       

kd = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2, 'X': 0, 'U': 0, 'B': 0, 'Z': 0}

#GRAVY: Grand average of hydropathicity index indicates the solubility of the proteins: 
#positive GRAVY (hydrophobic), negative GRAVY (hydrophilic) (Kyte and Doolittle, 1982).

# calculate the gravy according to kyte and doolittle.
def Gravy(ProteinSequence):

	if ProteinSequence.islower():
		ProteinSequence = ProteinSequence.upper()

	ProtGravy=0.0
	for i in ProteinSequence:
		ProtGravy += kd[i]

	return ProtGravy/len(ProteinSequence)
	
#----suppression numeros accession et full name si plusieurs-----

def calcul_length(child):
    end = 0
    begin = 0
    for location in child :
        for be in location :
            if be.tag == 'begin' :
                begin = float(be.attrib['position'])
            if be.tag == 'end' :
                end = float(be.attrib['position'])
        L = end - begin
    return L
    
for parent in tree.findall('.//entry'):
    taille = 0
    helixL = 0
    strandL = 0
    turnL = 0
    nb_accession=0
    nb_name=0
    nb_strand = 0   
    nb_turn = 0 
    nb_helix = 0
    tissue=[]
    cpt_cys=0
    for child in parent :
        if child.tag=='accession':
            nb_accession+=1
            if nb_accession>1:
                child.getparent().remove(child)
        if child.tag=='fullName':
            nb_name+=1
            if nb_name>1:
                child.getparent().remove(child) 
        if child.tag=='feature' and child.attrib['type']=='strand':
            nb_strand+=1
            strandL+=calcul_length(child)
        if child.tag=='feature' and child.attrib['type']=='turn':
            nb_turn+=1
            turnL+=calcul_length(child)
        if child.tag=='feature' and child.attrib['type']=='helix':
            nb_helix+=1
            helixL+=calcul_length(child)
            
        if child.tag=='tissue':
            if child.text not in tissue :
                tissue.append(child.text)
            else:
                child.getparent().remove(child)
        if child.tag =='sequence':
            sequence = child.text
            taille = float(child.attrib['length'])
            for aa in sequence :
                if aa=='C':
                    cpt_cys+=1
            sequence = sequence.replace('\n', '')
            gravy = Gravy(sequence)
            params  = ProtParam.ProteinAnalysis(sequence)
            phi = params.isoelectric_point()
            child.attrib["gravy"] = str(round(gravy,2))
            child.attrib["cysteine"] = str(cpt_cys)
            child.attrib["phi"] = str(round(phi,2))
      
    #print taille
    strandPct = (strandL/taille)*100
    helixPct = (helixL/taille)*100
    turnPct = (turnL/taille)*100
    strand = etree.Element("strand")
    strand.attrib["nb"] = str(nb_strand)
    strand.attrib["pct"] = str(round(strandPct,2))
    parent.append(strand)
    helix = etree.Element("helix")
    helix.attrib["nb"] = str(nb_helix)
    helix.attrib["pct"] = str(round(helixPct,2))
    parent.append(helix)
    turn = etree.Element("turn")
    turn.attrib["nb"] = str(nb_turn)
    turn.attrib["pct"] = str(round(turnPct,2))
    parent.append(turn)
    
    
    

tree.write("homocalc.xml")

