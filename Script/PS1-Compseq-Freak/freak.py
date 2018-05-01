# Freak
import sys
import random
import argparse
from Fonction import uniq,motif,sequence

def getparser():
	''' Cette fonction permet de prendre en compte les différentes options du programme'''
	parser = argparse.ArgumentParser(description='Calculate the composition of unique words in sequences')
	parser.add_argument('-seqall', action="store", dest='seq',type=str, required=True, help='Sequence(s) filename and optional format')
	parser.add_argument('-window', action="store", dest='W',default = 30, type=int, help='[30] Averaging window (Any integer value)')
	parser.add_argument('-step', action="store", dest='step',default = 1, type=int, help="[1] Stepping value (Any integer value)")
	parser.add_argument('-sreverse1', action="store", dest='sr1', type=bool,default = False,  help="Reverse (if DNA)")
	parser.add_argument('-sbegin1', action="store", dest='begin',type=int, default = 1, help='[1]Start of each sequence to be used')
	parser.add_argument('-send1', action="store", dest='end',type=int, default = 1, help='[1]End of each sequence to be used')
	parser.add_argument('-scircular1', action="store", dest='ci',type=bool,default = False, help='Y]Sequence is circular')
	parser.add_argument('-letters', action="store", dest='lt',type=str,required=True, help='[gc] Residue letters (Any string)')
	parser.add_argument('-stdout', action="store", dest='sd',type=bool,default = False, help='[Y]Write first file to standard output')
	return parser

	
def main() :
	''' Fonction qui permet, dans une séquence, de calculer les fréquences des résidus données (letters) dans un cadre de lecture donnée (window) avec un pas donné(step)'''
	parser=getparser()
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	Arguments= parser.parse_args()  # Permet pour la suite de donner une variable à l'argument pris en compte par la fonction command()
	fichier = Arguments.seq
	n = Arguments.W
	step = Arguments.step
	sreverse1 = Arguments.sr1
	sbegin1 =Arguments.begin
	send1 = Arguments.end
	circular = Arguments.ci
	letters = Arguments.lt
	letters = letters.lower()
	stdout = Arguments.sd
	reverse = False
	
	tableau = '' # On crée un tableau vide pour le compléter dans la suite du programme
	
	for i in sequence(fichier,n,reverse,sbegin1,send1,circular,sreverse1): # Permet de traiter les séquences une à une
		id,seq,seqRev = i # Permet de nommer les variables du tuple fait par la fonction sequence()
		if 'S' and 'R' and 'L' in seq : # Permet d'enlever les options qui ne sont pas necessaires pour des séquences protéiques
			sreverse1 = False
			circular = False
		seq = seq.lower() # Permet de normaliser le format de la séquence en mettant tous les résidus en minuscule
		seqRev = seqRev.lower()
		compteur = sbegin1 - step # Le compteur commence au premier nucléotide donné en prenant en compte le pas
		if send1 == 1 : # Permet de modifier l'entête du fichier si send1 n'est pas activée
			tableau = tableau+ 'FREAK of ' + id + ' from ' + str(sbegin1) + ' to ' + str(len(seq))+' Window ' +str(n)+ ' Step ' +str(step) +'\n' 
		else :  # Permet de modifier l'entete du fichier si send1 est activée
			tableau = tableau+ 'FREAK of ' + id + ' from ' + str(sbegin1) + ' to ' + str(send1) +' Window ' +str(n)+ ' Step ' +str(step) +'\n' 
		result = {} # Permet de créer un dictionnaire vide, que le programme utilisera plus tard pour le calcul de fréquence
		lletters = [] # Permet de créer une liste vide que le programme utilisera plus tard pour lister les résidus donnés en arguments
		for i in range(0,len(seq)-n+1,step): # Permet d'avoir une variable i qui permettra de créer le cadre de lecture et qui prend en compte le pas (step)
			frequence = 0 # Initialise la valeur de la fréquence recherchée à 0
			for lettre in letters : # Prend une par une chaque lettre (lettre) de la variable letters (tapée en arguments)
				lletters.append(lettre) # On ajoute une par une chaque lettre de l'option letters dans la liste
			lletters = uniq(lletters) # Permet d'enlever les doublons si l'utilisateur rentre des doublons dans l'option letters
			for j in lletters : # Ici la boucle permet de prendre chaque lettre de la liste faite précédemment
				result[j] = 0 # Crée dans le dictionnaire un couple dont la clé est la lettre et la valeur est 0
				for elt in seq[i:i+n]: # Permet de prendre chaque résidu de la séquence du cadre de lecture (cadre de lecture = seqRev[i:i+n])
					if elt == j: # Permet de sélectionner les éléments identiques à la lettre donnée
						result[j]+=1 # Permet d'incrémenter la valeur lorsque la condition précédente est vraie
				
				frequence = result[j]/n + frequence # Permet de calculer la fréquence de chaque lettre présente dans l'option letters dans le cadre de lecture
			
			compteur = compteur + step
			tableau =  tableau + '\n' + '%-10d %8.6f' % (compteur,frequence)
			result = {}
		fout = id +'.composition'
		tableau = tableau + '\n'
	if stdout == False :
		f = open(fout,'w')
		f.write(tableau)
		f.close()
	else :
		print(tableau)

main()
