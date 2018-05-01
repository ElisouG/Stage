import sys
import random

def uniq(liste):
	''' Permet de supprimé tous les doublons dans une liste '''
	liste = sorted(liste)
	liste2 = []
	liste2.append(liste[0])
	for i in range(len(liste)-1):
		if liste[i] != liste[i+1]:
			liste2.append(liste[i+1])
	return liste2

def motif(n,seq,ignorebz) :
	''' crée tous le motif possible de taille n '''
	adnl = []
	if 'S' and 'R' and 'L' in seq : # Permet de déterminer que la séquence est une séquence protéique
		if ignorebz == 'False' : # Si l'option ignorebz est fause alors on prend en compte aussi les acide aminés B et Z d'ou l'ajout dans notre liste
			base = 'ABCDEFGHIKLMNPQRSTUVWYZ'
		else :
			base='ACDEFGHIKLMNPQRSTUVWY'
	else : 
		base = 'ACGT'
	for i in base : # Suite de boucle qui permet de créer tous les motif dans l'ordre alphabétique puis que la variable base est rangé dans l'ordre alphabétique
		if n == 1 :
			adnl.append(i)
		else :
			for j in base:
				if n == 2 :
					adnl.append(i+j)
				else :
					for k in base :
						if n == 3 :
							adnl.append(i+j+k)
						else :
							for l in base :
								if n == 4 :
									adnl.append(i+j+k+l)
								else :
									for m in base :
										if n == 5 :
											adnl.append(i+j+k+l+m)
										else :
											for o in base :
												if n >= 6 :
													adnl.append(i+j+k+l+m+o)
													
	return adnl # retourne une liste de tous les motif possible

def seqembl(fichier):
    '''Permet d'extraire la séquence d'un fichier fasta'''
    ouvrir = open(fichier,'r')
    gene = ""
    seq = ""
    while True: # Boucle qui permet d'extraire un titre et une séquence d'un fichier ensembl
        ligne = ouvrir.readline()
        if ligne[0:2] == "AC":
            gene = ligne[2:]
            gene = gene.strip()
            gene = gene.replace(gene[-1],'')
        if ligne[0] == ' ':
            seq = seq + ligne[2:]
        if ligne[0:2] =="//":
            break
    ouvrir.close()
    chiffres = ( '0','1','2','3','4','5','6','7','8','9')
    for uneLettre in seq: # Permet de sélectionner que l'ID dans le tire de la séquence
        seq = seq.replace(" ","")
        seq = seq.replace("\n","")
        if uneLettre in chiffres:
            seq = seq.replace(uneLettre,'')
    seq = seq.upper()
    return [(gene, seq)]


def seqfasta(fichier):
	''' permet d'extraire la/les séquence(s) d'un fichier fasta'''
	fichier = fichier.strip()
	fichier = open(fichier,"r")
	resultat = []
	toutesLignes = fichier.readlines()
	fichier.close()
	x= 0
	seq=""
	for uneLigne in toutesLignes[0:] : # Boucle qui permet de crée une liste de tuple contenant l'id et la séquence
		if uneLigne[0] == '>' :
			if x == 0:
				ide = uneLigne.rstrip("\n").strip(">")
				ide = ide.replace(' ','|').split('|') # Permet de remplacer tous les espaces par un '|' pour le traitement de l'ID
				if '.' in ide[0] or '_' in ide[0] : # Permet de différencier si le fichier vients de NCBI ou ebi
					ide = ide[0] # Permet de récupérer l'id dans in fichier de NCBI
				else :
					ide = ide[2] # Permet de récupérer l'id dans in fichier de ebi
			else:
				resultat.append((ide,seq)) # ajout le tuple ide et seq dans notre liste
				ide = uneLigne.rstrip("\n").strip(">")
				ide = ide.replace(' ','|').split('|')
				if '.' in ide[0] or '_' in ide[0] :
					ide = ide[0]
				else :
					ide = ide[2]
				seq=""
		elif uneLigne != '\n':
			x+=1
			seq = seq + uneLigne.strip("\n")
	resultat.append((ide,seq))
	return resultat


def Reverse(seq) :
	''' Permet de faire le brin compélementaire de la séquence double brin'''
	seq = seq.upper()
	seqRev = seq.replace('A','t').replace('T','a').replace('C','g').replace('G','c')
	seqRev = seqRev.upper()
	seqRev1 = ''
	for i in range(1,len(seqRev)+1):
		seqRev1 =  seqRev1 + seqRev[-i] 
	return seqRev1



def sequence(fichier,n,reverse,sbegin1,send1,circular,sreverse1): 
	''' Fonction qui permet de traité la séquence en fonction de son format 
	et de la modifié en fonction des demande ( reverse, adn circulaire...)'''
	if 'fasta' in fichier : # Condition qui permet de traiter les fichier fasta
		sequences = seqfasta(fichier)
	else :
		sequences = seqembl(fichier)
	Resultat = []
	sbegin1 = sbegin1 -1
	for i in sequences: # Cette boucle permet de traiter chaque tuple de la liste (id, séquence) une par une
			id,seq = i
			seq = seq.upper()
			if sreverse1 == True :
				seq = Reverse(seq)
			if circular == True :
				seq = seq + seq[0:n-1] # Permet de rajouter les résidus nécéssaires (pour le motif ou le cadre de lecture) du début de la séquence à la fin de la séquence
			if send1 != 1 :
				seq = seq[0:send1] # Permet de finir la séquence au résidu indiqué par l'option send1
			seq = seq.replace(seq[:sbegin1],'',sbegin1) # Permet de commencer au résidu indiqué par l'option sbegin1
			seqRev = ''
			if reverse == True :
				seqRev = Reverse(seq)# Permet de faire la séquence complémentaire de notre séquence
			Resultat.append((id,seq,seqRev)) # Créer une liste de tuple à 3 éléments contenant l'id, la séquence et la séquence complémentaire	
	return Resultat
	
