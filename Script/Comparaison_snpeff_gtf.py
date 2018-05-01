#!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Comparaison_snpeff_gtf.py
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

	pathGTF = "/media/sf_DATA/Stage_UM-ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"

	GTF = open(pathGTF) # Lire le fichier GTF initial
	linesGTF = GTF.readlines() # Sauvegarder dans la variable linesGTF les lignes du GTF initial
	GTF.close() # Fermer le fichier
	compa = open("/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/GTF_vs_snpEff.txt", "w") # Création d'un fichier
	compa.write('%20s\t%20s\t%20s\t%20s\n' %("nbCDS","nbexon","nbstart","nbstop")) # Ecrit les entêtes du fichier
	nbCDS = 0 # Initialise le compteur CDS à 0
	nbstop = 0 # Initialise le compteur Stop à 0
	nbexon = 0 # Initialise le compteur exon à 0
	nbstart = 0 # Initialise le compteur start à 0
	for line in linesGTF: # Boucle qui va permettre de compter les occurences des compteurs précédemment initialisé
		nbCDS = nbCDS + line.count('CDS') # Compte le nombre de CDS présent
		nbexon = nbCDS + line.count('exon') # Compte le nombre d'exon présent
		nbstart = nbCDS + line.count('start_codon') # Compte le nombre de start présent
		nbstop = nbCDS + line.count('stop_codon') # Compte le nombre de stop présent
	compa.write('%20s\t%20s\t%20s\t%20s\n' %(nbCDS,nbexon,nbstart,nbstop)) # Ecrit dans le tableau les valeurs obtenues pour chacun des compteurs
	compa.close() # Ferme le fichier
