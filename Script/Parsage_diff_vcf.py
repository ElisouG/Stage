#!/home/egueret/miniconda3/bin/python3
##!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Parsage_diff_vcf.py
# @author Elise GUERET


if __name__ == "__main__":
	# Fichiers existants
	pathAll = "/media/sf_DATA/filtering/Compa_filtered.snps_All.diff.sites_in_files"
	pathC = "/media/sf_DATA/filtering/Compa_filtered.snps_C.diff.sites_in_files"
	pathJ = "/media/sf_DATA/filtering/Compa_filtered.snps_J.diff.sites_in_files"
	# Fichiers créés
	pathF1 = "/media/sf_DATA/filtering/Common_Sites_All.txt"
	pathF2 = "/media/sf_DATA/filtering/Common_Sites_C.txt"
	pathF3 = "/media/sf_DATA/filtering/Common_Sites_J.txt"

	################# Fichier All ####################

	all_sites = open(pathAll,"r") 
	allLines = all_sites.readlines() 
	all_sites.close()

	f1 = open(pathF1,"w") # Création d'un nouveau fichier
	f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  | brin \n' %('chromosome','old_start','old_end','new_start','new_end','codon','Length add')) # Nomme les colonnes du fichier
	f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))
	

	for line in allLines:
		if "B" in line:
			print(line)



	################# Fichier C ####################
	
	c_sites = open(pathC,"r") 
	cLines = c_sites.readlines() 
	c_sites.close()

	f2 = open(pathF2,"w") # Création d'un nouveau fichier
	f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  | brin \n' %('chromosome','old_start','old_end','new_start','new_end','codon','Length add')) # Nomme les colonnes du fichier
	f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))
	

	################# Fichier J ####################
	
	j_sites = open(pathJ,"r") 
	jLines = j_sites.readlines() 
	j_sites.close()

	f3 = open(pathF3,"w") # Création d'un nouveau fichier
	f3.write('%15s | %15s | %15s | %15s | %15s | %15s  | brin \n' %('chromosome','start','end','first_frame','last_frame','brin')) # Nomme les colonnes du fichier
	f3.write('%15s | %15s | %15s | %15s | %15s | %15s  | ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))