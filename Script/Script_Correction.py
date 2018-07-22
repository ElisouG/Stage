#!/home/egueret/miniconda3/bin/python3
##!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Correction_Annotation.py
# @author Elise GUERET


#################### Commande lancement sur le cluster ###########################

#qsub -cwd -V -S /bin/bash -l h_rt=99:00:00 -M elise.gueret@gmail.com -m bes -N Test_sge_modification_gtf -o /home/egueret/Stage_UM_ISEM/sge_test_modification_gtf.out -e /home/egueret/Stage_UM_ISEM/sge_test_modification_gtf.err -b y "python3 /home/egueret/Stage/Script/Correction_annotation.py"

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
from Fonction_Correction import *


if __name__ == "__main__":

#################### Path File Ordi bureau   ######################

	# pathFasta = "/media/sf_DATA/Stage_UM_ISEM/Puce_57K/genome_IGV/labrax.fa"
	# pathWarning = "/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err"
	# pathGTF = "/media/sf_DATA/Stage_UM_ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"
	# pathAnnotation = "/media/sf_DATA/Stage_UM-ISEM/Puce_57K/correction_annotation.gtf"
	# pathF2 = "/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/NoStop_exo.txt"
	# pathF1 = "/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/NoStart_exo.txt"
	# pathAnntotationNew = "/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/New_Annotation_20_04_2018.gtf"
	# PathGenomeR = ''
	# pathTranscriptome = "/media/sf_DATA/Stage_UM_ISEM/Puce_57K/Dicentrarchus_merged-transcript.gtf"
	# pathTranscrits = "/media/sf_DATA/Stage_UM_ISEM/SnpEff_last/Liste_transcrits.txt"

#################### Path File Ordi_portable   ######################

	# pathFasta = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax.fa"
	# pathWarning = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err"
	# pathGTF = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"
	# pathAnnotation = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/correction_annotation.gtf"
	# pathF2 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/NoStop_exo.txt"
	# pathF1 = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/NoStart_exo.txt"
	# pathAnntotationNew = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/SnpEff_last/New_Annotation_20_04_2018.gtf"
	# PathGenomeR = '/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/genome_IGV/labrax_reverse.fa'
	
#################### Path File Cluster    ######################

	pathFasta = "/home/egueret/Stage_UM_ISEM/Puce_57K/genome_IGV/labrax.fa"
	pathWarning = "/home/egueret/Stage_UM_ISEM/SnpEff_last/sge_test_snpEff_puce_last_debug.err"
	pathGTF = "/home/egueret/Stage_UM_ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"
	pathAnnotation = "/home/egueret/Stage_UM_ISEM/Puce_57K/correction_annotation.gtf"
	pathF2 = "/home/egueret/Stage_UM_ISEM/SnpEff_last/NoStop_exo_18_05_18.txt"
	pathF1 = "/home/egueret/Stage_UM_ISEM/SnpEff_last/NoStart_exo_18_05_18.txt"
	pathF3 = "/home/egueret/Stage_UM_ISEM/SnpEff_last/MultipleStop_13_06_18.txt"
	pathF4 = "/home/egueret/Stage_UM_ISEM/SnpEff_last/Incomplete_13_06_18.txt"
	pathAnntotationNew = "/home/egueret/Stage_UM_ISEM/SnpEff_last/New_Annotation_18_05_2018.gtf"
	PathGenomeR = '/home/egueret/Stage_UM_ISEM/Puce_57K/genome_IGV/labrax_reverse.fa'
	pathTranscriptome = "/home/egueret/Stage_UM_ISEM/Puce_57K/Dicentrarchus_merged-transcript.gtf"
	pathTranscrits = "/home/egueret/Stage_UM_ISEM/SnpEff_last/Liste_transcrits.txt"



	pathStopTested = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/NoStop_Tested_13_06_18.txt"
	pathStartTested ="/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/NoStart_Tested_13_06_18.txt"
	pathCDS = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/CDS_13_06_18.txt"
	pathSequenceCDS = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/Sequence_CDS_13_06_18.txt"
	pathProteinesCDS = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/Protéines_CDS_13_06_18.txt"

################### Récupération info #################

	listeInfo = waringParse(pathWarning)

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

	listeNoStart,listeNoStop = correctionSens(listeInfo,GenomeR,Genome,pathF1,pathF2,pathF3,pathF4)		
				
	#################### Comparaison snpEff et GTF ########################
	print('Création fichier NoStart/NoStop')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))

########################Pasage transcriptome ######################

	listeTranscrit = transcriptomeParse(pathTranscriptome)

###################### Modification du GTF #############################

	print('Debut Annotation')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))


	newAnnotation(pathAnntotationNew,pathAnnotation,listeNoStop)
	
	print('Fin Annotation')
	print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))

######################## Vérification dans le trancriptome #########################
	
	listeNoStopTested = TestStopCorrection(pathStopTested,listeNoStop,listeTranscrit)
	listeNoStartTested = TestStopCorrection(pathStartTested,listeNoStart,listeTranscrit)

######################## Recherche des positions des CDS #########################
	
	CDSFinaux = recupPosCDS(pathAnntotationNew,pathCDS)

	######################## Recherche des séquences des CDS #########################
	
	CDSComplete = recupSeqCDS(pathSequenceCDS,CDSFinaux)



	######################## TRADUCTION DES CDS EN PROTÉINES #########################

	# SequenceCDS = open(pathSequenceCDS)
	# linesSequenceCDS= SequenceCDS.readlines()
	# SequenceCDS.close()

	ProteinesCDS = open(pathProteinesCDS, "w")
	ProteinesCDS.write("%s | %s | %s | %s\n" % ('Chromosome','geneID','Filtre','Protéines'))

	for elt in CDSComplete:
		listeErreur = []
		K1 = elt[0]
		geneID = elt[1]
		seqNt = elt[3]
		if seqNT[0:2] == 'ATG' and seqNT[-3:] in ['TAA','TAG','TGA'] and len(seqNT)%3 == 0 :
			if K1 != 'MT':
				seqProt = str(seqNt.translate())
				nbreStar = seqProt.count('*')
				if nbreStar > 1:
					Filtre = 'Plusieurs stop' 
					ProteinesCDS.write("%s | %s | %s | %s\n" % (K1,geneID,Filtre,seqProt))	
				elif nnbreStar == 1:
					Filtre = 'RAS'
					ProteinesCDS.write("%s | %s | %s | %s\n" % (K1,geneID,Filtre,seqProt))
				elif nbreStar == 0:
					Filtre = 'Pas de stop'
					ProteinesCDS.write("%s | %s | %s | %s\n" % (K1,geneID,Filtre,seqProt))
			elif K1 == 'MT':
				seqProt = str(seqNt.translate(table="Vertebrate Mitochondrial",))
				nbreStar = seqProt.count('*')
				if nbreStar > 1:
					Filtre = 'Plusieurs stop' 
					ProteinesCDS.write("%s | %s | %s | %s\n" % (K1,geneID,Filtre,seqProt))	
				elif nbreStar == 1:
					Filtre = 'RAS'
					ProteinesCDS.write("%s | %s | %s | %s\n" % (K1,geneID,Filtre,seqProt))
				elif nbreStar == 0:
					Filtre = 'Pas de stop'
					ProteinesCDS.write("%s | %s | %s | %s\n" % (K1,geneID,Filtre,seqProt))	
		else :
			erreur = []
			if  seqNT[0:2] != 'ATG' :
				erreur.append('Start_erreur')
			
			if  seqNT[-3:] not in ['TAA','TAG','TGA'] :
				erreur.append('Stop_erreur')
			if  slen(seqNT)%3 != 0  :
				erreur.append('Len_erreur')

			listeErreur.append(elt.append(erreur))
	ProteinesCDS.close()

	for elt in listeErreur :
		print(elt)