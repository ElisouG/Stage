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
	######################## Vérification des Nostart et Nostop ######################

	# Lecture du fichier transcriptome et sauvegarde des lignes dans linesTranscriptome
	Transcriptome = open(pathTranscriptome)
	linesTranscriptome = Transcriptome.readlines()
	Transcriptome.close()

	listeTranscrit = [] # Création d'une liste pour y stocker les transcrits
	# Parsage du transcriptome
	Transcrits = open(pathTranscrits, "w")
	for line in linesTranscriptome:
		K = line.split('\t')[0]
		S = line.split('\t')[3]
		E = line.split('\t')[4]
		tID = line.split('\t')[8].split('"')[3]
		if line.split('\t')[8].split('"')[3] != tID:
			listeTranscrit.append([K,S,E,tID])
		elif line.split('\t')[8].split('"')[3] == tID:
			E = line.split('\t')[4]
		Transcrits.write(listeTranscrit)
	Transcrits.close()
