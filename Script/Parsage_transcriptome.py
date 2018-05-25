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
	a = "none"
	listeNoStopTested = []
	for elt in listeNoStop:
		chr2 = elt[0]
		st2 = elt[1]
		en2 = elt[2]
		for elt in listeTranscrit:
			K = elt[0]
			S = elt[1]
			E = elt[2]
			if chr2 == K: # Vérifier égalité des chromosomes
				if S > st2 > E and S > en2 < E:
					a = "PASS"
				elif S < st2 > E and S > en2 > E::
					a = "Not Valid"
	listeNoStopTested.append([chr2,st2,en2,a])

	b = "none"
	listeNoStartTested = []
	for elt[] in listeNoStart:
		chr1 = elt[0]
		st1 = elt[1]
		en1 = elt[2]
		for elt[] in listeTranscrit:
			K = elt[0]
			S = elt[1]
			E = elt[2]
			if chr1 == K: # Vérifier égalité des chromosomes
				if S > st1 > E and S > en1 < E:
					b = "PASS"
				elif S < st1 > E and S > en1 > E::
					b = "Not Valid"
	listeNoStopTested.append([chr1,st1,en1,b])





	######################## TRADUCTION DES CDS EN PROTÉINES #########################

	prot = mySeq.translate() # http://www.python-simple.com/python-biopython/manipulation-sequences.php 