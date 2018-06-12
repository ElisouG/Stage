#!/home/egueret/miniconda3/bin/python3
##!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Parsage_transcriptome.py
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
	pathTranscrits = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/Liste_transcrits.txt"

	pathF2 = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_17-05-18/NoStop_exo_18_05_18.txt"
	pathF1 = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_17-05-18/NoStart_exo_18_05_18.txt"

	pathStopTested = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/NoStop_Tested_12_06_18.txt"
	pathStartTested ="/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/NoStart_Tested_12_06_18.txt"

	pathAnntotationNew = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_17-05-18/New_Annotation_18_05_2018.gtf"
	pathCDS = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/CDS_12_06_18.txt"

	pathFasta = "/home/egueret/Stage_UM_ISEM/Puce_57K/genome_IGV/labrax.fa"
	PathGenomeR = '/home/egueret/Stage_UM_ISEM/Puce_57K/genome_IGV/labrax_reverse.fa'
	pathSequenceCDS = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/Sequence_CDS_12_06_18.txt"
	pathProteinesCDS = "/home/egueret/Stage_UM_ISEM/SnpEff_Modif_28-05-18/Protéines_CDS_12_06_18.txt"

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
				Transcrits.write("%s\t%s\t%s\t%s\n" % (K,S,E,tID))
				listeTranscrit.append([K,S,E,tID])
				
			S = line.split('\t')[3]
			K = line.split('\t')[0]
			E = line.split('\t')[4]
		elif line.split('\t')[8].split('"')[3] == tID:
			E = line.split('\t')[4]
		tID = line.split('\t')[8].split('"')[3]
	listeTranscrit.append([K,S,E,tID])
	Transcrits.write("%s\t%s\t%s\t%s\n" % (K,S,E,tID))
	Transcrits.close()

	######################## Parsage NoStart ######################

	NoStart = open(pathF1)
	linesNoStart = NoStart.readlines()
	NoStart.close()
	linesNoStart = linesNoStart[2:len(linesNoStart)]

	listeNoStart = []


	for line in linesNoStart:
		chr1 = line.replace(' ','').split('|')[0]
		st1 = line.replace(' ','').split('|')[3]
		en1 = line.replace(' ','').split('|')[4]
		listeNoStart.append([chr1,st1,en1])

	######################## Parsage NoStop ######################
	
	NoStop = open(pathF2)
	linesNoStop = NoStop.readlines()
	NoStop.close()
	linesNoStop = linesNoStop[2:len(linesNoStop)]

	listeNoStop = []


	for line in linesNoStop:
		chr2 = line.replace(' ','').split('|')[0]
		st2 = line.replace(' ','').split('|')[3]
		en2 = line.replace(' ','').split('|')[4]
		listeNoStop.append([chr2,st2,en2])

	######################## Vérification dans le trancriptome #########################
	
	StopTested = open(pathStopTested, "w")
	StopTested.write('%s\t%s\t%s\t%s\t%s\n' %('chromosome','CDS_start','CDS_end','new_end','Filter'))
	listeNoStopTested = []

	for elt in listeNoStop:
		Pass = False
		chr2 = elt[0]
		st2 = int(elt[1])
		en2 = int(elt[2])
		for elt in listeTranscrit:
			K = elt[0]
			S = int(elt[1])
			E = int(elt[2])
			if chr2 == K: # Vérifier égalité des chromosomes
				if S < en2 <= E: # and S < st2 < E:
					listeNoStopTested.append([chr2,st2,en2,"PASS"])
					StopTested.write("%s\t%s\t%s\t%s\t%s\n" % (chr2,S,E,en2,"PASS"))
					Pass = True
					break
		if Pass == False :
			listeNoStopTested.append([chr2,st2,en2,"Not Valid"])
			StopTested.write("%s\t%s\t%s\t%s\t%s\n" % (chr2,S,E,en2,"Not Valid"))
	StopTested.close()
	

	StartTested = open(pathStartTested, "w")
	StartTested.write('%s\t%s\t%s\t%s\t%s\n' %('chromosome','CDS_start','CDS_end','new_start','Filter'))
	listeNoStartTested = []
	for elt in listeNoStart:
		Pass = False
		chr1 = elt[0]
		st1 = int(elt[1])
		en1 = int(elt[2])
		for elt in listeTranscrit:
			K = elt[0]
			S = int(elt[1])
			E = int(elt[2])
			if chr1 == K: # Vérifier égalité des chromosomes
				if S <= st1 < E: # and S < en1 < E:
					listeNoStartTested.append([chr1,st1,en1,"PASS"])
					StartTested.write("%s\t%s\t%s\t%s\t%s\n" % (chr1,S,E,st1,"PASS"))
					Pass = True
					break
		if Pass == False :
			listeNoStartTested.append([chr1,st1,en1,"Not Valid"])
			StartTested.write("%s\t%s\t%s\t%s\t%s\n" % (chr1,S,E,en1,"Not Valid"))
	StartTested.close()


	######################## Recherche des positions des CDS #########################
	
	NewGTF = open(pathAnntotationNew)
	linesNewGTF= NewGTF.readlines()
	NewGTF.close()

	listeCDS = []
	CDSFinaux = []

	CDS = open(pathCDS, "w")
	CDS.write("%s | %s | %s | %s | %s:%s\n" % ('Chromosome','geneID','brin','frame','CDS_start','CDS_end'))
	geneID = "none"
	for line in linesNewGTF:
		if 'CDS' in line and geneID == "none":
			K1 = line.split('\t')[0]
			brin = line.split('\t')[6]
			CDS_start = line.split('\t')[3]
			CDS_end = line.split('\t')[4]
			frame = line.split('\t')[7]
			geneID = line.replace('"','').split()[1]
			CDS.write("%s | %s | %s | %s | %s:%s" % (K1,geneID,brin,frame,CDS_start,CDS_end)) 
			listeCDS = [K1,geneID,brin,frame,'%s:%s'% (CDS_start,CDS_end)]
		elif line.replace('"','').split()[1] != geneID :
			geneID = line.replace('"','').split()[1]
			K1 = line.split('\t')[0]
			brin = line.split('\t')[6]
			CDS_start = line.split('\t')[3]
			CDS_end = line.split('\t')[4]
			frame = line.split('\t')[7]
			CDS.write("%s | %s | %s | %s | %s:%s\n" % (K1,geneID,brin,frame,CDS_start,CDS_end)) 
			listeCDS = [K1,geneID,brin,frame,'%s:%s'% (CDS_start,CDS_end)] 
			CDSFinaux.append(listeCDS)
		elif line.replace('"','').split(' ')[1] == geneID :
			CDS_start = line.split('\t')[3]
			CDS_end = line.split('\t')[4]
			CDS.write(" | %s:%s " % (CDS_start,CDS_end))
			listeCDS.append('%s:%s'% (CDS_start,CDS_end))  
			geneID = line.replace('"','').split(' ')[1]
	CDS.close()

	######################## Recherche des séquences des CDS #########################

	Genome = fasta2dict(pathFasta)
	GenomeR = fasta2dict(PathGenomeR)

	#SequenceCDS = open(pathSequenceCDS, "w")
	#SequenceCDS.write("%s | %s | %s | %s\n" % ('Chromosome','geneID','brin','Sequence'))

	CDSComplete = []
	for elt in CDSFinaux:
		seqFinale = ""
		print(elt)
		K1 = elt[0]
		geneID = elt[1]
		brin = elt[2]
		frame = elt[3]
		positionCDS=CDSFinaux[4:len(CDSFinaux)]
		for elt in positionCDS: # Est-ce qu'il ne faudrait pas que j'échange cette ligne avec celle du dessous?
			if brin == '+':
				sequence = Genome[chromosome].seq 
				if frame == '1':
					CDS_start = str(int((elt.split(':')[0]))+1)
					CDS_end = elt.split(':')[1]
					seqCDS = sequence[CDS_start:CDS_end] 
					seqFinale = seqFinale+seqCDS
					#SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
				elif frame == '2':
					CDS_start = str(int((elt.split(':')[0]))+2)
					CDS_end = elt.split(':')[1]
					seqCDS = sequence[CDS_start:CDS_end]
					seqFinale = seqFinale+seqCDS
					#SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
				elif frame == '0':
					CDS_start = elt.split(':')[0]
					CDS_end = elt.split(':')[1]
					seqCDS = sequence[CDS_start:CDS_end]
					seqFinale = seqFinale+seqCDS
					#SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
			
			elif brin == '-':
				sequence = GenomeR[chromosome].seq	
				if frame == '1':
					CDS_start = str(int((elt.split(':')[0]))+1)
					CDS_end = elt.split(':')[1]
					seqCDS = sequence[CDS_start:CDS_end]
					seqFinale = seqFinale+seqCDS
					#SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
				elif frame == '2':
					CDS_start = str(int((elt.split(':')[0]))+2)
					CDS_end = elt.split(':')[1]
					seqCDS = sequence[CDS_start:CDS_end]
					seqFinale = seqFinale+seqCDS
					#SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
				elif frame == '0':
					CDS_start = elt.split(':')[0]
					CDS_end = elt.split(':')[1]
					seqCDS = sequence[CDS_start:CDS_end]
					seqFinale = seqFinale+seqCDS
					#SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
					seqFinale = seqFinale+seqCDS
	CDSComplete.append(K1,geneID,brin,seqFinale)	
	#SequenceCDS.close()

	######################## TRADUCTION DES CDS EN PROTÉINES #########################

	# SequenceCDS = open(pathSequenceCDS)
	# linesSequenceCDS= SequenceCDS.readlines()
	# SequenceCDS.close()

	ProteinesCDS = open(pathProteinesCDS, "w")
	ProteinesCDS.write("%s | %s | %s | %s\n" % ('Chromosome','geneID','Filtre','Protéines'))

	for elt in CDSComplete:
		K1 = elt[0]
		geneID = elt[1]
		seqNt = elt[4]
		if K1 != 'MT':
			seqProt = str(seqNt.translate(cds=True))
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
			seqProt = str(seqNt.translate(table="Vertebrate Mitochondrial",cds=True))
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
	ProteinesCDS.close()
