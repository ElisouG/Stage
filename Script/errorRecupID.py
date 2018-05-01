#!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package errorRecupID.py
# @author Elise GUERET

if __name__ == "__main__":
	file = open("/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err","r")# Ouvrir le fichier à lire
	lines = file.readlines()# Mettre dans lines toutes les lignes du fichier ci-dessus
	file.close()# Fermer le fichier pour ne pas qu'il reste ouvert
	f = open("/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/warningRecupID.txt","w")# Création d'un nouveau fichier
	f.close()# Fermer le fichier créé pour ne pas qu'il reste ouvert
	a = 0# Initialisation de a.
	for line in lines:# Boucle qui va permettre de récupérer le message d'erreur ainsi que la position sur le génome
		if line[0:7] == "WARNING" or line[0:5] == "ERROR":# Pour trouver les débuts de lignes
			line = line.split("'")[-1]
			a = 1
			message = "\tMessage d'erreur : " + line
		elif a == 1 : 
			a = 0
			line = line.split("|")[0]
			f = open("/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/warningRecupID.txt","a")
			f.write(line + message+'\n')
			f.close()
