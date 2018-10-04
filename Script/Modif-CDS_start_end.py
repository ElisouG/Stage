#!/home/egueret/miniconda3/bin/python3
##!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Modif-CDS_start_end.py
# @author Elise GUERET


#qsub -cwd -V -S /bin/bash -l h_rt=99:00:00 -M elise.gueret@gmail.com -m bes -N Test_sge_Script_Correction -o /home/egueret/Stage_UM_ISEM/sge_test_Script_Correction.out -e /home/egueret/Stage_UM_ISEM/sge_test_Script_Correction.err -b y "python3 /home/egueret/Stage/Script/Modif-CDS_start_end.py"

if __name__ == "__main__":

#################### Path File Cluster    ######################

	pathAnnotation = "/home/egueret/Stage_UM_ISEM/Puce_57K/correction_annotation.gtf"

################### Récupération info #################

	info = open(pathAnnotation,"r")
	lines = info.readlines() # Sauvegarder dans la variable lines les lignes du fichier précédent
	info.close() # Fermer le fichier

	listeFct = [] # Création d'une liste vide
	listeFinale = [] # Création d'une liste vide
	a = False # Initialisation du compteur
	
	for line in lines:
		lineSplit = line.split('\t')
		K = lineSplit[0]
		brin = lineSplit[6]
		gene_id = lineSplit[8].split(" ")[1].replace('"',"")
		fct = lineSplit[2]
		st = lineSplit[3]
		end = lineSplit[4]
		listeFct.append([gene_id,brin,K,fct,st,end])

		for elt in listeFct :
			if gene_id == lineSplit[8].split(" ")[1].replace('"',"") :
				gene_id = gene_id
				brin = brin
				k = K 	
				if brin == '+' and fct == "CDS" or fct == "exon" and a = False :
					a = True
					fct = lineSplit[2]
					st = lineSplit[3]
					end = lineSplit[4]
				elif brin == '-' and fct == "CDS" or fct == "exon" and a = True:
					a = False
					fct = lineSplit[2]
					st = lineSplit[3]
					end = lineSplit[4]
		listeFinale.append([gene_id,brin,K,fct,st,end])
	print(listeFinale)