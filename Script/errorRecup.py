#!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package errorRecup.py
# @author Elise GUERET

if __name__ == "__main__":

    #################### Initialisation variable ###################
	nonLast = False 
	stop=False
	start=False
	mutiple=False # Initialisation de a

	#################### Path File ######################
	pathF1 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/warningStop.txt"
	pathF2 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/warningStart.txt"
	pathF3 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/warningincomplete.txt"
	pathF4 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/warningMultipleStop.txt"

	###################### open file #####################
	file = open("/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err","r")# Ouvrir le fichier à lire
	lines = file.readlines()# Mettre dans lines toutes les lignes du fichier ci-dessus
	file.close()# Fermer le fichier pour ne pas qu'il reste ouvert
	print(lines[-1])
	# f1 = open(pathF1,"w")# Création d'un nouveau fichier
	# f1.close()# Fermer le fichier créé pour ne pas qu'il reste ouvert
	# f2 = open(pathF2,"w")# Création d'un nouveau fichier
	# f2.close()# Fermer le fichier créé pour ne pas qu'il reste ouvert
	# f3 = open(pathF3,"w")# Création d'un nouveau fichier
	# f3.close()# Fermer le fichier créé pour ne pas qu'il reste ouvert
	# f4 = open(pathF4,"w")# Création d'un nouveau fichier
	# f4.close()# Fermer le fichier créé pour ne pas qu'il reste ouvert


	# ##################### Main ##################################
	# for line in lines:# Boucle qui va permettre de récupérer le message d'erreur ainsi que la position sur le génome
	# 	if line[0:7] == "WARNING" or line[0:5] == "ERROR":# Pour trouver les débuts de lignes
	# 		type = line.split("'")[-1]
	# 		if "non-last" in type:
	# 			textNl= "\tMessage d'erreur : " + type
	# 			nonLast= True
	# 		elif "STOP" in type:
	# 			textStop="\tMessage d'erreur : " + type
	# 			stop = True
	# 		elif "START" in type:
	# 			textStart="\tMessage d'erreur : " + type
	# 			start =True
	# 		elif "mutiple" in type:
	# 			textMutiple="\tMessage d'erreur : " + type
	# 			mutiple = True
	# 	else:
	# 		if nonLast == True : 
	# 			info = line.split("|")[0]
	# 			f4 = open(pathF4,"a")
	# 			f4.write(info + textNl)
	# 			f4.close()
	# 			nonLast=False
	# 		if stop : 
	# 			info = line.split("|")[0]
	# 			f1 = open(pathF1,"a")
	# 			f1.write(info + textStop)
	# 			f1.close()
	# 			stop=False
	# 		if start : 
	# 			info = line.split("|")[0]
	# 			f2 = open(pathF2,"a")
	# 			f2.write(info + textStart)
	# 			f2.close()
	# 			start=False
	# 		if mutiple : 
	# 			info = line.split("|")[0]
	# 			f3 = open(pathF3,"a")
	# 			f3.write(info + textMutiple)
	# 			f3.close()
	# 			mutiple=False