#!/home/egueret/miniconda3/bin/python3
##!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Correction_Annotation.py
# @author Elise GUERET


#################### Commande lancement sur le cluster ###########################

# qsub -cwd -V -S /bin/bash -l h_rt=99:00:00 -M elise.gueret@gmail.com -m bes -N Test_sge_parsage -o /home/egueret/Stage_UM_ISEM/sge_test_parsage.out -e /home/egueret/Stage_UM_ISEM/sge_test_parsage.err -b y "python3 /home/egueret/Stage/Script/Parsage_transcriptome.py"

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
	pathTranscrits = "/home/egueret/Stage_UM_ISEM/SnpEff_last/Liste_transcrits.txt"

	pathF2 = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_17-05-18/NoStop_exo_18_05_18.txt"
	pathF1 = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_17-05-18/NoStart_exo_18_05_18.txt"

	######################## Parsage Transcriptome ######################

	# Lecture du fichier transcriptome et sauvegarde des lignes dans linesTranscriptome
	Transcriptome = open(pathTranscriptome)
	linesTranscriptome = Transcriptome.readlines()
	Transcriptome.close()

	listeTranscrit = [] # Création d'une liste pour y stocker les transcrits
	
	Transcrits = open(pathTranscrits, "w")
	tID = "none"
	for line in linesTranscriptome:
		if line.split('\t')[8].split('"')[3] != tID :
			if tID != "none" :
				listeTranscrit.append([K,S,E,tID])
				Transcrits.write("%s\t%s\t%s\t%s" % (K,S,E,tID))
			S = line.split('\t')[3]
			K = line.split('\t')[0]
			E = line.split('\t')[4]
		elif line.split('\t')[8].split('"')[3] == tID:
			E = line.split('\t')[4]
		tID = line.split('\t')[8].split('"')[3]
	listeTranscrit.append([K,S,E,tID])
	Transcrits.write("%s\t%s\t%s\t%s" % (K,S,E,tID))
	Transcrits.close()

	######################## Parsage NoStart ######################

	NoStart = open(pathF1)
	linesNoStart = NoStop.readlines()
	NoStop.close()

	listeNoStart = []

	chr1 = "none"
	st1 = "none"
	en1 = "none"
	for line in linesNoStop:
		chr1 = line.split('|')[0]
		st1 = line.split('|')[3]
		en1 = line.split('|')[4]
	listeNoStart.append([chr1,st1,en1])

	######################## Parsage NoStop ######################
	
	NoStop = open(pathF2)
	linesNoStop = NoStop.readlines()
	NoStop.close()

	listeNoStop = []

	chr2 = "none"
	st2 = "none"
	en2 = "none"
	for line in linesNoStop:
		chr2 = line.split('|')[0]
		st2 = line.split('|')[3]
		en2 = line.split('|')[4]
	listeNoStop.append([chr2,st2,en2])

	######################## Vérification dans le trancriptome #########################

	for elt in listeNoStop:
		for elt in listeTranscrit:
			if elt[0] in listeNoStop == elt[0] in listeTranscrit: # Vérifier égalité des chromosomes
				if elt[1] in listeTranscrit > elt[1] in listeNoStop > elt[2] in listeTranscrit and elt[1] in listeTranscrit > elt[2] in listeNoStop > elt[2] in listeTranscrit:


	for elt[] in listeNoStart:
		for elt[] in listeTranscrit:
			if elt[0] in listeNoStop == elt[0] in listeTranscrit: # Vérifier égalité des chromosomes
				if elt[1] in listeTranscrit > elt[1] in listeNoStop > elt[2] in listeTranscrit and elt[1] in listeTranscrit > elt[2] in listeNoStop > elt[2] in listeTranscrit:




	######################## TRADUCTION DES CDS EN PROTÉINES #########################

	prot = mySeq.translate() # http://www.python-simple.com/python-biopython/manipulation-sequences.php 