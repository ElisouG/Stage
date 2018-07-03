#!/home/egueret/miniconda3/bin/python3
##!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Correction_multiples.py
# @author Elise GUERET


#################### Commande lancement sur le cluster ###########################

# qsub -cwd -V -S /bin/bash -l h_rt=99:00:00 -M elise.gueret@gmail.com -m bes -N Test_sge_parsage -o /home/egueret/Stage_UM_ISEM/sge_test_parsage.out -e /home/egueret/Stage_UM_ISEM/sge_test_parsage.err -b y "python3 /home/egueret/Stage/Script/Correction_multiples.py"

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


if __name__ == "__main__":

	pathTranscriptome = "/home/egueret/Stage_UM_ISEM/Puce_57K/Dicentrarchus_merged-transcript.gtf"
	pathGTF = "/home/egueret/Stage_UM_ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"

	######################## Récupération frame GTF ###########################

	GTF = open(pathGTF)
	linesGTF = GTF.readlines()
	GTF.close()

	ListeExons = []

	for line in linesGTF:
		if 'exon' in line:
			K = line.split('\t')[0]
			S = line.split('\t')[3]
			E = line.split('\t')[4]
			F = line.split('\t')[7]
			ListeExons.append([K,S,E,F])
			#print(ListeExons)		
		




	######################## Récupération frame transcriptome ###########################

	Transcriptome = open(pathTranscriptome)
	linesTranscriptome = Transcriptome.readlines()
	Transcriptome.close()

	ListeTranscrits = []

	for line in linesTranscriptome:
		K = line.split('\t')[0]
		S = line.split('\t')[3]
		E = line.split('\t')[4]
		F = line.split('\t')[7]
		ListeTranscrits.append([K,S,E,F])
		#print(ListeTranscrits)
