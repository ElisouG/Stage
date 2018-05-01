#!/usr/bin/python3
## !/root/miniconda3/bin/python3
#-*- coding: utf-8 -*-
# @package Correction_Annotation.py
# @author Elise GUERET


##################### Importation Module #######################
## Python modules
import argparse, os, subprocess, sys, time, glob, re
from time import localtime, strftime, sleep, clock, time

## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet, IUPAC, Gapped
from Bio import AlignIO
from Bio.Nexus import Nexus

##### module perso #######
from fonctions import seqfasta, Reverse
from module_seb import fasta2dict

if __name__ == "__main__":

	#################### Path File ######################

	# pathFasta = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax.fa"
	# pathWarning = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err"
	# pathGTF = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"
	# pathAnnotation = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/correction_annotation.gtf"

	# pathF2 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/NoStop.txt"
	# pathF1 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/NoStart.txt"
	pathFasta = "/media/sf_DATA/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax.fa"
	pathWarning = "/media/sf_DATA/Stage_UM-ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err"
	pathGTF = "/media/sf_DATA/Stage_UM-ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"
	pathAnnotation = "/media/sf_DATA/Stage_UM-ISEM/Puce_57K/correction_annotation.gtf"

	pathF2 = "/media/sf_DATA/Stage_UM-ISEM/SnpEff_last/NoStop.txt"
	pathF1 = "/media/sf_DATA/Stage_UM-ISEM/SnpEff_last/NoStart.txt"	

	################### Récupération info #################

	info = open(pathWarning,"r") # Lire le fichier warning obtenu par snpEff
	lines = info.readlines() # Sauvegarder dans la variable lines les ligne du fichier précédent
	info.close() # Fermer le fichier
	a = 0 # Initialisation du compteur
	b = 0 # Initialisation du compteur
	pb = [] # Création d'une liste vide
	listeInfo = [] # Création d'une seconde liste vide

	for line in lines:
		# Récupérer le type de problème
		if line[0:8] == "WARNING:" or line[0:6] == "ERROR:" : # Permet de récupérer le début de la ligne
			# Triage des problèmes en fonction de leur origine
			a = 1 # Mettre le compteur a à 1
			if "non-last" in line.split("'")[-1] : # Récupérer les multiple stop
				pb.append("Multiple Stop") # Ajouter ce pb à la liste des pb
			elif "STOP" in line.split("'")[-1] : # Récupérer les no stop
				pb.append("No Stop") # Ajouter ce pb à la liste des pb
			elif "START" in line.split("'")[-1] : # Récupérer les no start
				pb.append("No Start") # Ajouter ce pb à la liste des pb
			elif "has length" in line.split("'")[-1] : # Récupérer les incomplets
				pb.append("Incomplet") # Ajouter ce pb à la liste des pb
		# Récupérer leur cadre de lecture
		elif "rank:" in line and b!=1:
			frame1 = line.split(',')[2].split(':')[1] # Récupérer le 1er cadre de lecture
			b=1 # Mettre le compteur b à 1
		elif b == 1 and "rank:" in line :
			frame2 = line.split(',')[2].split(':')[1] # Récupérer le second cadre de lecture
		# Récupérer le reste des informations
		elif a == 1 and len(pb) == 1: # Vérifier si le compteur est à 1
			a = 0 # Ré-initialisation du compteur
			lineSplit = line.split(",") # Séparation des caractères de la ligne par ","
			if "-" in lineSplit[1] : # Récupération de l'information brin anti-sens
				brin = "anti-sens"
			elif "+" in lineSplit[1] : # Récupération de l'information brin sens
				brin = "sens"
			lineSplit = lineSplit[0].split(":") # Séparation des caractères de la ligne par ":"
			chromosome = lineSplit[0] # Récupération du chromosome
			lineSplit = lineSplit[1].split("-") # Séparation des caractères de la ligne par "-"
			start = int(lineSplit[0]) # Récupération du début du gène et mise en integer (chiffre)
			end = int(lineSplit[1]) # Récupération de la fin du gène et mise en integer (chiffre)
		elif line[0:5] == '\t\tCDS' and len(pb) == 1:
			b = 0 # Ré-initialisation du compteur
			for elt in pb :
				listeInfo.append([chromosome,start,end,brin,elt,frame1,frame2]) # Création d'une liste comportant des listes des différentes info sur les problèmes rencontrés
			pb = [] # Ré-initialisation de la liste
		elif a == 1 and len(pb) != 1:
			pb = [] # Ré-initialisation de la liste

	# for elt in listeInfo:
	# 	print(elt)
	Genome = fasta2dict(pathFasta) # Crée un dictionnaire du génome

	############ Creation séquence reverse #########################
	#listeId = Genome.keys() # Récupérer le nom des chromosomes sous forme d'une liste
	#f = open('/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax_reverse.fa','w') # Création du fichier du génome reverse
	#for elt in listeId:
	#	f.write('>%s\n%s\n'%(elt,Genome[elt].seq.reverse_complement())) # Ecrit comme dans un fichier fasta le chromosome puis la séquence
	#f.close() # Fermer le fichier
	#GenomeR = fasta2dict('/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax_reverse.fa') # Crée un dictionnaire du génome reverse précédemment créé
	GenomeR = fasta2dict('/media/sf_DATA/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax_reverse.fa') # Crée un dictionnaire du génome reverse précédemment créé

	################### Recherche codons start, stop ##################

	f1 = open(pathF1,"w") # Création d'un nouveau fichier
	f1.write('%20s\t%20s\t%20s\t%20s\t%20s\n' %('chromosome','gene_start','début_codon','fin_codon','codon')) # Nomme les colonnes du fichier
	f2 = open(pathF2,"w") # Création d'un nouveau fichier
	f2.write('%20s\t%20s\t%20s\t%20s\t%20s\n' %('chromosome','gene_end','début_codon','fin_codon','codon')) # Nomme les colonnes du fichier
	codonStop = ['TAA','TAG','TGA'] # Liste des différents codons stop existants
	codonStart = 'ATG' # Nomme le codon start
	codonStart2 = 'AAG' # Nomme le codon start alternatif
	listeNoStop = []
	for elt in listeInfo :
		# Récupération des différents éléments de la liste listeInfo
		chromosome = elt[0] # Récupération du chromosome
		start = elt[1] # Récupération de la position du début du gène
		start1 = start - 100 # Ajout de 100bp en amont du gène
		end = elt[2] # Récupération de la position de la fin du gène
		end1 = end + 999 # Ajout de 1000bp en aval du gène
		brin = elt[3] # Récupération du brin
		pb = elt[4] # Récupération du problème
		if brin == 'sens':
			sequence = Genome[chromosome].seq # Si le brin est sens alors il prend le fichier "sens"
			if pb == 'No Stop' :
				for i in range(end+1,end1-2,3): # Recherche des codons stop en aval du gène
					if sequence[i:i+3] in codonStop :
						f2 = open(pathF2,"a") # Ouverture du fichier "No Stop"
						f2.write('%20s\t%20s\t%20s\t%20s\t%20s\n' %(chromosome,end-1,i+1,i+3,sequence[i:i+3])) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
						f2.close() # Fermer le fichier
						listeNoStop.append([chromosome,str(int(end)-1),str(int(end)+1),str(i+1),str(i+3)])
						break # S'arrête dès qu'il trouve le premier codon stop


		elif brin == 'anti-sens':
			sequence = GenomeR[chromosome].seq # Si le brin est anti-sens alors il prend le fichier "anti-sens"
			length = len(sequence)
			if pb == 'No Stop' :
				for i in range(length-(start+1),length-(start -1002),3): # Recherche des codons stop en aval du gène
					if sequence[i:i+3] in codonStop :
						f2 = open(pathF2,"a") # Ouverture du fichier "No Stop"
						f2.write('%20s\t%20s\t%20s\t%20s\t%20s\n' %(chromosome,start,str(length-(i+2)+1),str(length-i+1),sequence[i:i+3])) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
						f2.close() # Fermer le fichier
						listeNoStop.append([chromosome,str(int(start)+1),str(int(start)+3),str(length-(i+2)+1),str(length-i+1)])
						break # S'arrête dès qu'il trouve le premier codon stop

		elif pb == 'No Start' :
			st = 'NA'
			en = 'NA'
			for i in range(start1,start-2,3): # Recherche des codons start en amont du gène
				if sequence[i:i+3] == codonStart :
					st = i
					en = i+2
			if st != 'NA':
				f1 = open(pathF1,"a") # Ouverture du fichier "No Start"
				f1.write('%20s\t%20s\t%20s\t%20s\t%20s\n' %(chromosome,start,st,en,sequence[st:en+1])) # Ecrire dans ce fichier les codons start trouvés ainsi que leur position
				f1.close() # Fermer le fichier 
				
	#################### Comparaison snpEff et GTF ########################

	# GTF = open(pathGTF) # Lire le fichier GTF initial
	# linesGTF = GTF.readlines() # Sauvegarder dans la variable linesGTF les lignes du GTF initial
	# GTF.close() # Fermer le fichier
	# compa = open("/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/GTF_vs_snpEff.txt", "w") # Création d'un fichier
	# compa.write('%20s\t%20s\t%20s\t%20s\n' %("nbCDS","nbexon","nbstart","nbstop")) # Ecrit les entêtes du fichier
	# nbCDS = 0 # Initialise le compteur CDS à 0
	# nbstop = 0 # Initialise le compteur Stop à 0
	# nbexon = 0 # Initialise le compteur exon à 0
	# nbstart = 0 # Initialise le compteur start à 0
	# for line in linesGTF: # Boucle qui va permettre de compter les occurences des compteurs précédemment initialisé
	# 	nbCDS = nbCDS + line.count('CDS') # Compte le nombre de CDS présent
	# 	nbexon = nbCDS + line.count('exon') # Compte le nombre d'exon présent
	# 	nbstart = nbCDS + line.count('start_codon') # Compte le nombre de start présent
	# 	nbstop = nbCDS + line.count('stop_codon') # Compte le nombre de stop présent
	# compa.write('%20s\t%20s\t%20s\t%20s\n' %(nbCDS,nbexon,nbstart,nbstop)) # Ecrit dans le tableau les valeurs obtenues pour chacun des compteurs
	# compa.close() # Ferme le fichier


	###################### Modification du GTF #############################
	
	Annotation = open(pathAnnotation)
	linesAnnotation = Annotation.readlines()
	Annotation.close()

	newAnnotation = open("/media/sf_DATA/Stage_UM-ISEM/SnpEff_last/New_Annotation.gtf", "w")
	# newAnnotation = open("/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/New_Annotation.gtf", "w")
	for line in linesAnnotation:
		for elt in listeNoStop:
			if elt[0] in line and '\t'+elt[1]+'\t' in line and 'stop_codon' in line: 
				line = line.replace(elt[1],elt[3])
				line = line.replace(elt[2],elt[4])
		newAnnotation.write(line)
	newAnnotation.close()




	# Prendre le chromosome, les positions du gènes dans le fichier .err ainsi que le brin sens ou antisens
	# Il faut un fichier fasta sens et un fichier fasta antisens
	# Chercher la séquence à ces positions dans le fichier fasta sens ou antisens
	# Extraire cette séquence avec 100 caractères supplémentaires 5'UTR et 1000 en 3'UTR
	# Savoir si cette séquence contient des N si cette séquence contient des N alors il faut éliminer la séquence
	# Grâce au cadre de lecture (CDS) chercher les codons start "ATG" ou alternatif "AAG" ou les codons Stop "TAA", "TAG", "TGA"


###################### Modification du GTF #############################	
	# Liste1 = ['LG4','1402956','1402957','1402959','TGA']
	# Liste2 = ['LG12','3449403','3449506','3449508','TAA']
	# Liste3 = ['LG19','4944694','4944785','4944787','TAA']
	# Liste4 = ['LG4','15640582','15640610','15640612','TAA']
	# Liste5 = ['UN','28485802','28485803','28485805','TGA']
	# Liste6 = ['LG16','17487889','17487917','17487919','TGA']
	# ListeNoStop = [Liste1,Liste2,Liste3,Liste4,Liste5,Liste6]
	# for elt in ListeNoStop:
	# 	elt[1] = str(int(elt[1])+1)
	# 	elt[2] = str(int(elt[2])+1)
	# 	elt[3] = str(int(elt[3])+1)
	# NoStop = open(pathF2)
	# linesNoStop = NoStop.readlines()
	# NoStop.close()