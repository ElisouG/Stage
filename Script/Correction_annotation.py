#!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Correction_Annotation.py
# @author Elise GUERET


#################### Commande lancement sur le cluster ###########################

# qsub -cwd -V -S /bin/bash -M elise.gueret@gmail.com -m bes -N Test_sge_test_modification_gtf -o /home/egueret/Stage_UM_ISEM/sge_test_modification_gtf.out -e /home/egueret/Stage_UM_ISEM/sge_test_modification_gtf.err -b y "python3 /home/egueret/Stage_UM_ISEM/script/Correction_annotation.py"

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

#################### Path File Ordi bureau   ######################

	# pathFasta = "/media/sf_DATA/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax.fa"
	# pathWarning = "/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err"
	# pathGTF = "/media/sf_DATA/Stage_UM_ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"
	# pathAnnotation = "/media/sf_DATA/Stage_UM-ISEM/Puce_57K/correction_annotation.gtf"
	# pathF2 = "/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/NoStop_exo.txt"
	# pathF1 = "/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/NoStart_exo.txt"
	# pathAnntotationNew = "/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/New_Annotation_20_04_2018.gtf"
	# PathGenomeR = ''

	#################### Path File Ordi_portable   ######################

	pathFasta = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax.fa"
	pathWarning = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err"
	pathGTF = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"
	pathAnnotation = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/correction_annotation.gtf"
	pathF2 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/NoStop_exo.txt"
	pathF1 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/NoStart_exo.txt"
	pathAnntotationNew = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/New_Annotation_20_04_2018.gtf"
	PathGenomeR = '/mnt/c/Users/missl/Documents//Stage_UM-ISEM/Puce_57K/genome_IGV/labrax_reverse.fa'
	
	
    #################### Path File Cluster    ######################

	# pathFasta = "/home/egueret/Stage_UM_ISEM/Puce_57K/genome_IGV/labrax.fa"
	# pathWarning = "/home/egueret/Stage_UM_ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err"
	# pathGTF = "/home/egueret/Stage_UM_ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"
	# pathAnnotation = "/home/egueret/Stage_UM_ISEM/Puce_57K/correction_annotation.gtf"
	# pathF2 = "/home/egueret/Stage_UM_ISEM/SnpEff_last/NoStop_exo_30_04_18.txt"
	# pathF1 = "/home/egueret/Stage_UM_ISEM/SnpEff_last/NoStart_exo_30_04_18.txt"
	# pathAnntotationNew = "/home/egueret/Stage_UM_ISEM/SnpEff_last/New_Annotation_30_04_2018.gtf"
	# PathGenomeR = '/home/egueret/Stage_UM_ISEM/Puce_57K/genome_IGV/labrax_reverse.fa'


	################### Récupération info #################

	info = open(pathWarning,"r") # Lire le fichier warning obtenu par snpEff
	lines = info.readlines() # Sauvegarder dans la variable lines les lignes du fichier précédent
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
			frame2 = line.split(',')[2].split(':')[1] # Récupérer le dernier cadre de lecture
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
	# 		print(elt)
	print('Warning Parse')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	Genome = fasta2dict(pathFasta) # Crée un dictionnaire du génome
	print('Genome')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))

	############ Creation séquence reverse #########################
	#listeId = Genome.keys() # Récupérer le nom des chromosomes sous forme d'une liste
	#f = open('/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax_reverse.fa','w') # Création du fichier du génome reverse
	#for elt in listeId:
	#	f.write('>%s\n%s\n'%(elt,Genome[elt].seq.reverse_complement())) # Ecrit comme dans un fichier fasta le chromosome puis la séquence
	#f.close() # Fermer le fichier
	GenomeR = fasta2dict(PathGenomeR) # Crée un dictionnaire du génome reverse précédemment créé
	
	print('GenomeR')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	################### Recherche codons start, stop ##################

	f1 = open(pathF1,"w") # Création d'un nouveau fichier
	f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  | brin \n' %('chromosome','old_start','old_end','new_start','new_end','codon','Length add')) # Nomme les colonnes du fichier
	f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))
	f2 = open(pathF2,"w") # Création d'un nouveau fichier
	f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  | brin \n' %('chromosome','old_start','old_end','new_start','new_end','codon','Length add')) # Nomme les colonnes du fichier
	f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))
	codonStop = ['TAA','TAG','TGA'] # Liste des différents codons stop existants
	codonStart = 'ATG' # Nomme le codon start
	codonStart2 = 'AAG' # Nomme le codon start alternatif
	listeNoStop = []
	listeNoStart = []
	for elt in listeInfo :
		# Récupération des différents éléments de la liste listeInfo
		chromosome = elt[0] # Récupération du chromosome
		startGTF = elt[1]
		start = elt[1] # Récupération de la position du début du gène
		start1 = start - 99 # Ajout de 100bp en amont du gène
		endGTF = elt[2]
		end = elt[2] # Récupération de la position de la fin du gène
		end1 = end + 999 # Ajout de 1000bp en aval du gène
		brin = elt[3] # Récupération du brin
		pb = elt[4] # Récupération du problème
		frame1 = elt[5] # Récupération du 1er cadre de lecture
		frame2 = elt[6] # Récupération du dernier cadre de lecture


#################################################Sens #######################################

		if brin == 'sens':
			sequence = Genome[chromosome].seq # Si le brin est sens alors il prend le fichier "sens"
			if pb == 'No Stop' :
				if frame2 == '1' :
					end = end - 2
					end1 = end1 - 2
				elif frame2 == '2' :
					end = end - 1 
					end1 = end1 -1
				for i in range(end+1,end1-2,3): # Recherche des codons stop en aval du gène
					if sequence[i:i+3] in codonStop :
						Longueur = (i+3) - (int(endGTF)+1)
						f2 = open(pathF2,"a") # Ouverture du fichier "No Stop"
						f2.write('%15s | %15s | %15s | %15s | %15s | %15s | %15s | +\n' %(chromosome,endGTF-1,endGTF+1,i+1,i+3,sequence[i:i+3],str(Longueur))) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
						f2.close() # Fermer le fichier
						listeNoStop.append([chromosome,str(int(endGTF)-1),str(int(endGTF)+1),str(i+1),str(i+3),'+'])
						break # S'arrête dès qu'il trouve le premier codon stop

			elif pb == 'No Start':
				st = 'NA'
				en = 'NA'
				if frame1 == '1':
					start = start +1
					start1 = start1 +1
				elif frame1 == '2':
					start = start + 2
					start1 = start1 + 2

				for i in range(start1,start-3,3): # Recherche des codons start en amont du gène
					if sequence[i:i+3] == codonStart :
						st = i
						en = i+2
				if st != 'NA':
					Longueur = (int(startGTF)+1) - int((st+1))
					f1 = open(pathF1,"a") # Ouverture du fichier "No Start"
					f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  +\n' %(chromosome,startGTF+1,startGTF+3,st+1,en+1,sequence[st:en+1],str(Longueur))) # Ecrire dans ce fichier les codons start trouvés ainsi que leur position
					f1.close() # Fermer le fichier 	
					listeNoStart.append([chromosome,str(int(startGTF)+1),str(int(startGTF)+3),str(st+1),str(en+1),'+'])
		

############################# Anti sens ##############################################
		elif brin == 'anti-sens':
			sequence = GenomeR[chromosome].seq # Si le brin est anti-sens alors il prend le fichier "anti-sens"
			length = len(sequence)
			if pb == 'No Stop' :
				if frame1 == '1' :
					start = start + 2
				elif frame1 == '2' :
					start = start +1
				for i in range(length-1-(start-1),length-1-(start-1002),3): # Recherche des codons stop en aval du gène
					if sequence[i:i+3] in codonStop :
						Longueur = (int(startGTF)+1) - (length-(i+1))
						f2 = open(pathF2,"a") # Ouverture du fichier "No Stop"
						f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  -\n' %(chromosome,startGTF +1,startGTF+3,str(length-(i+3)),str(length-(i+1)),sequence[i:i+3],str(Longueur))) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
						f2.close() # Fermer le fichier
						listeNoStop.append([chromosome,str(int(startGTF)+1),str(int(startGTF)+3),str(length-(i+3)),str(length-(i+1)),'-'])
						break # S'arrête dès qu'il trouve le premier codon stop

			elif pb == 'No Start' :
				st = 'NA'
				en = 'NA'
				if frame2 == '1':
					end = end -1
				if frame2 == '2':
					end = end -2
				for i in range(length-1-(end+99),length-1-(end +1) ,3): # Recherche des codons start en amont du gène
					if sequence[i:i+3] == codonStart :
						st = i
						en = i+2
				if st != 'NA':
					Longueur = (length-(st)) - (int(endGTF)+1)
					f1 = open(pathF1,"a") # Ouverture du fichier "No Start"
					f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  -\n' %(chromosome,endGTF-1,endGTF+1,length-(en),length-(st),sequence[st:en+1],str(Longueur))) # Ecrire dans ce fichier les codons start trouvés ainsi que leur position
					f1.close() # Fermer le fichier 
					listeNoStart.append([chromosome,str(int(endGTF)-1),str(int(endGTF)+1),str(length-(en)),str(length-(st)),'-'])
				
	#################### Comparaison snpEff et GTF ########################
	print('Création fichier NoStart/NoStop')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))


	###################### Modification du GTF #############################
	print('Debut Annotation')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	Annotation = open(pathAnnotation)
	linesAnnotation = Annotation.readlines()
	Annotation.close()

	print('ListeAnnotation Done')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	#newAnnotation = open("/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/New_Annotation_20_04_2018.gtf", "w")
	newAnnotation = open(pathAnntotationNew, "w")
	for line in linesAnnotation:
		for elt in listeNoStop:
			if elt[0] in line and '\t'+elt[1]+'\t' in line and 'stop_codon' in line: 
				line = line.replace(elt[1],elt[3])
				line = line.replace(elt[2],elt[4])
			if  elt[5] == '+'  and  elt[0] in line:
				if 'exon' in line and elt[2] in line  :
					line = line.replace(elt[2],elt[4])
				elif 'CDS' in line and str(int(elt[2])-3) in line :
					line = line.replace(str(int(elt[2])-3),str(int(elt[4])-3))
			if  elt[5] == '-' and  elt[0] in line : 
				if 'exon' in line and elt[1] in line :
					line = line.replace(elt[1],elt[3])
				elif 'CDS' in line and str(int(elt[1])+3) in line:
					line = line.replace(str(int(elt[1])+3),str(int(elt[3])+3))
		for elt in listeNoStart:
			if elt[0] in line and '\t'+elt[1]+'\t' in line and 'start_codon' in line: 
				line = line.replace(elt[1],elt[3])
				line = line.replace(elt[2],elt[4])
			if  elt[5] == '+'  and  elt[0] in line:
				if elt[1] in line  :
					line = line.replace(elt[1],elt[3])
			if  elt[5] == '-' and  elt[0] in line: 
				if elt[2] in line :
					line = line.replace(elt[2],elt[4])

		newAnnotation.write(line)
	newAnnotation.close()
	print('Fin Annotation')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))

