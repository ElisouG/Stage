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

		
	Transcrits = open(pathTranscrits, "w")
	tID = "none"
	for line in linesTranscriptome:
		if tID == "none" :
			tID = line.split('\t')[8].split('"')[3]
			S = line.split('\t')[3]
			E = line.split('\t')[4]
			Transcrits.write("%15s | %15s:%15s" % (tID,S,E))
		elif line.split('\t')[8].split('"')[3] != tID :
			tID = line.split('\t')[8].split('"')[3]
			S = line.split('\t')[3]
			E = line.split('\t')[4]
			Transcrits.write("\n %15s | %15s:%15s" % (tID,S,E))
		elif line.split('\t')[8].split('"')[3] == tID:
			S = line.split('\t')[3]
			E = line.split('\t')[4]
			Transcrits.write("| %15s:%15s " % (S,E))
		tID = line.split('\t')[8].split('"')[3]
	Transcrits.close()