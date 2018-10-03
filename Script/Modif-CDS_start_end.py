#!/home/egueret/miniconda3/bin/python3
##!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Modif-CDS_start_end.py
# @author Elise GUERET


qsub -cwd -V -S /bin/bash -l h_rt=99:00:00 -M elise.gueret@gmail.com -m bes -N Test_sge_Script_Correction -o /home/egueret/Stage_UM_ISEM/sge_test_Script_Correction.out -e /home/egueret/Stage_UM_ISEM/sge_test_Script_Correction.err -b y "python3 /home/egueret/Stage/Script/Modif-CDS_start_end.py"

if __name__ == "__main__":

#################### Path File Cluster    ######################

	pathAnnotation = "/home/egueret/Stage_UM_ISEM/Puce_57K/correction_annotation.gtf"

################### Récupération info #################

	info = open(pathAnnotation,"r")
	lines = info.readlines() # Sauvegarder dans la variable lines les lignes du fichier précédent
	info.close() # Fermer le fichier

	listeFct = [] # Création d'une seconde liste vide
	
	for line in lines:
		K = line.split("	")[0]
		fct = line.split("	")[2]
		st = line.split("	")[3]
		end = line.split("	")[4]
		brin = line.split("	")[6]
		gene_id = line.split("	")[8].split(" ")[1].replace('"',"")
		listeFct.append([gene_id,brin,K,fct,st,end])
	print(listeFct)