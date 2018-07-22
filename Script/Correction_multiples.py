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
			B = line.split('\t')[6]
			F = line.split('\t')[7]
			gID = line.split('\t')[8].split('"')[1].split('|')[0]
			ListeExons.append([K,S,E,B,F,gID])
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
		B = line.split('\t')[6]
		F = line.split('\t')[7]
		gID = line.split('\t')[8].split('"')[7]
		ListeTranscrits.append([K,S,E,B,F,gID])
		print(ListeTranscrits)


	###################### Egalité frame ? Recherche frame #########################

	ListegIDEgaux = []

	for elt in ListeExons:
		K = elt[1]
		S = elt[2]
		E = elt[3]
		B = elt[4]
		F = elt[5]
		gID = elt[6]
		for elt in ListeTranscrits:
			K1 = elt[1]
			S1 = elt[2]
			E1 = elt[3]
			B1 = elt[4]
			F1 = elt[5]
			gID1 = elt[6]
			if gID == gID1:
				T = 'égaux'
				ListegIDEgaux.append([gID,T])
				if F == F1 and F == '.':
					Tester tous les frames
				elif F == F1 and F == '1':
					Tester les frames 2 et 3
				elif F == F1 and F == '2':
					Tester les frames 1 et 3
				elif F == F1 and F == '3':
					Tester les frames 1 et 2
			elif gID != gID1:
				T = 'différents'
				ListegIDEgaux.append([gID,gID1,T])


				
