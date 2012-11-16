import lxml.etree as etree

file = open('echantillon_parse.xml','r')
tree = etree.parse(file)
root = tree.getroot()

#------recup mots cles--------

kw = ['Cell adhesion', 'Transcription regulation', 'Repeat', 'RNA-binding', 'Developmental protein', 'Early protein', 'Zinc-finger', 'Cell division', 'Viral envelope protein', 'Protease', 'Transcription', 'Differentiation', 'Phosphoprotein', 'Immunoglobulin domain', 'Cell cycle', 'Reference proteome', 'Complete proteome', 'Viral nucleoprotein', '3D-structure', 'Ubl conjugation', 'Mitochondrion', 'Viral attachment to host cell', 'Polymorphism', 'Lipid metabolism', 'Ion transport', 'Nucleotidyltransferase', 'Methylation', 'Transferase', 'Acetylation', 'Ribonucleoprotein', 'Palmitate', 'Activator', 'Viral penetration into host cytoplasm', 'Cell membrane', 'Coiled coil', 'Secreted', 'G-protein coupled receptor', 'Ubl conjugation pathway', 'Olfaction', 'AIDS', 'Signal', 'Immunity', 'Transport', 'Inhibition of host interferon signaling pathway by virus', 'Host cell membrane', 'Apoptosis', 'Fusion of virus membrane with host endosomal membrane', 'Transit peptide', 'Host membrane', 'Sensory transduction', 'Cytoplasm', 'Signal-anchor', 'Membrane', 'Fusion of virus membrane with host membrane', 'Cell junction', 'Oxidoreductase', 'Repressor', 'Virion', 'Late protein', 'Myristate', 'Cytoskeleton', 'Calcium', 'Glycoprotein', 'Serine/threonine-protein kinase', 'Cell projection', 'Viral immunoevasion', 'Capsid protein', 'Nucleotide-binding', 'Protein transport', 'mRNA processing', 'Magnesium', 'Ionic channel', 'Receptor', 'Metal-binding', 'Transmembrane', 'DNA replication', 'Host cytoplasm', 'ATP-binding', 'Disease mutation', 'Viral RNA replication', 'Golgi apparatus', 'Kinase', 'Isopeptide bond', 'Transmembrane helix', 'Transducer', 'Inhibition of host innate immune response by virus', 'Disulfide bond', 'Zinc', 'Nucleus', 'Alternative splicing', 'Direct protein sequencing', 'Lipoprotein', 'Host nucleus', 'Hydrolase', 'Cleavage on pair of basic residues', 'Virus entry into host cell', 'Host-virus interaction', 'DNA-binding', 'Endoplasmic reticulum']





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

# ------supprime les mots cles non selectionnes-----		
for node in tree.getiterator():
	if node.tag=='keyword':
		if node.text not in kw :
			#print node.text
			node.getparent().remove(node)
			

tree.write("Test.xml")

