#!/home/egueret/miniconda3/bin/python3
##!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Exo_Flo.py
# @author Elise GUERET


#################### Commande lancement sur le cluster ###########################

#qsub -cwd -V -S /bin/bash -l h_rt=99:00:00 -M elise.gueret@gmail.com -m bes -N Test_sge_parsage -o /home/egueret/Stage_UM_ISEM/sge_test_parsage.out -e /home/egueret/Stage_UM_ISEM/sge_test_parsage.err -b y "python3 /home/egueret/Stage/Script/Parsage_transcriptome.py"

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

	pathTranscriptome = "/home/egueret/Stage_UM_ISEM/Puce_57K/Dicentrarchus_merged-transcript.gtf"
	pathTranscrits = "/home/egueret/Stage_UM_ISEM/Liste_transcrits_Exo_Flo.txt"



	######################## Parsage Transcriptome ######################

	# Lecture du fichier transcriptome et sauvegarde des lignes dans linesTranscriptome
	Transcriptome = open(pathTranscriptome)
	linesTranscriptome = Transcriptome.readlines()
	Transcriptome.close()
	sousliste = []
	listeFinal = []
	Transcrits = open(pathTranscrits, "w")
	tID = "none"
	for line in linesTranscriptome:
		if tID == "none" :
			tID = line.split('\t')[8].split('"')[3]
			S = line.split('\t')[3]
			E = line.split('\t')[4]
			# Transcrits.write("%s | %s:%s" % (tID,S,E)) # Faux car tu vas avoir un ligne avec None a la place du nom de transcript
		elif line.split('\t')[8].split('"')[3] != tID :
			tID = line.split('\t')[8].split('"')[3]
			S = line.split('\t')[3]
			E = line.split('\t')[4]
			Transcrits.write("\n %s | %s:%s" % (tID,S,E))
			listeFinal.append(sousliste)
			sousliste = [tID,%s:%s% (S,E)] # Tu initialise la sous liste qui va juste prendre une fois le Tid 
					               # et a chaque fois ajouté les positions
		elif line.split('\t')[8].split('"')[3] == tID:
			S = line.split('\t')[3]
			E = line.split('\t')[4]
			Transcrits.write(" | %s:%s " % (S,E))
			sousliste.append(%s:%s% (S,E))  
		tID = line.split('\t')[8].split('"')[3]
	# Tu as oublié la dernière ligne comme (regarde le script que l'on a fait tous les deux, comme le 
	# script récupère le transcript d'avant le dernier est pas récupéré :p
	Transcrits.close()

	#Transcrits = open(pathTranscrits)
	#linesTranscrits = Transcrits.readlines()
	#Transcrits.close()

	#listeTranscrits = []
	#for line in linesTranscrits:
		#line = line.replace("|",",")
		#listeTranscrits.append([line])


