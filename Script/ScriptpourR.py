#!/usr/bin/python3.5
#-*- coding: utf-8 -*-
# @package ScriptpourR.py
# @author Elise GUERET


##################### Importation Module #######################
## Python modules
import argparse, os, subprocess, sys, time, glob, re
from time import localtime, strftime, sleep, clock, time


# RScript bqsr.R ${sampleName}_recalibration.csv ${sampleName}_marked_duplicates_sorted_recal_data_before.table ${sampleName}_plots_RAnalyzeCovariates.pdf

path = "/media/bigvol/egueret/Donnees_CRECHE/R_results/"


for nameFIle in os.listdir(path):
	if nameFIle.endswith('_recalibration.csv') : # La tu récupère juste les csv
		nameFIle = nameFile.replace('_recalibration.csv','') # La tu vas juste récupérer le {sampleName} que tu voulais si j'ai bien compris
		pathFile = path+nameFile (attention nameFile est que le nom du fichier et non tous le chemin^^)
		print('RScript bqsr.R %s_recalibration.csv  %s_marked_duplicates_sorted_recal_data_before.table %s_plots_RAnalyzeCovariates.pdf'%(pathFile,pathFile,pathFile))