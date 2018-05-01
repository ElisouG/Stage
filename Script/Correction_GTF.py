#!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Correction_GTF.py
# @author Elise GUERET


################## Importation des modules #######################

import os

################# Fonctions ##############

def extractWarning(file):  # Fonction rapide
	'''Cette fonction permet d'extraire les données des messages d'erreur données par SnpEff
	La fonction renvoie une liste contenant des sous-listes. Chaque sous-liste représente une erreur 
	et contient le chromosome, le début du gène, la fin du gène, le brin, le cadre de lecture du 1er exon, 
	le cadre de lecture du dernier exon et le message d'erreur.

	@file : Fichier d'erreurs de SnpEff
	'''
	f = open(file,"r")
	listeLines = f.readlines()
	f.close()
	a = False # Permet de récupérer le deuxième ligne d'erreur
	b = False # Permet de récupérer frame1
	c = False # Permet de récupérer frame2
	listePb = []
	listeInfo = []
	for line in listeLines:
		if "ERROR:" == line[0:6] or "WARNING:" == line[0:8]:
			message = line.split("'")[-1]
			if 'non-last' in message:
				listePb.append("Multiple Stop")
			elif 'not a start' in message:
				listePb.append("No start")
			elif 'not a stop' in message:
				listePb.append("No stop")
			elif 'has length' in message:
				listePb.append("Incomplete")
			a = True
		elif a:
			a = False 
			line = line.split(",")
			geneStart = line[0].split(":")[1].split("-")[0]
			geneEnd = line[0].split(":")[1].split("-")[1]
			strand = line[1].split(":")[1]
			chromosome = line[0].split(":")[0]
		elif "Exons:" in line : 
			b = True
		elif b == True and "rank:" in line :
			b = False
			frame1 = line.split("frame:")[1].split(",")[0]
			c = True
		elif c == True and "rank:" in line :
			frame2 = line.split("frame:")[1].split(",")[0]	
		elif "CDS" in line and len(listePb) == 1:
			listeInfo.append([chromosome, int(geneStart), int(geneEnd), strand, int(frame1), int(frame2), listePb[0]])
			listePb = []
			c = False
		elif "CDS" in line and len(listePb) >= 2:
			print(listePb)
			listePb = []
			c = False
	return listeInfo




			














if __name__ == '__main__':
	file = 'media/sf_DATA/Stage_UM_ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err'
	listeInfo = extractWarning(file)
	for elt in listeInfo:
		print(elt)