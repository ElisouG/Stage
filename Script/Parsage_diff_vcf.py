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
	pathF1B = "/media/sf_DATA/filtering/Exclude_Sites_All.txt"
	pathF2 = "/media/sf_DATA/filtering/Common_Sites_C.txt"
	pathF2B = "/media/sf_DATA/filtering/Exclude_Sites_C.txt"
	pathF3 = "/media/sf_DATA/filtering/Common_Sites_J.txt"
	pathF3B = "/media/sf_DATA/filtering/Exclude_Sites_J.txt"

	################# Fichier All ####################

	all_sites = open(pathAll,"r") 
	allLines = all_sites.readlines() 
	all_sites.close()

	f1 = open(pathF1,"w") # Création d'un nouveau fichier
	f1.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %('CHROM','POS1','POS2','IN_FILE','REF1','REF2','ALT1','ALT2')) # Nomme les colonnes du fichier
	f1.close()
	f1B = open(pathF1B,"w") # Création d'un nouveau fichier
	f1B.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %('CHROM','POS1','POS2','IN_FILE','REF1','REF2','ALT1','ALT2')) # Nomme les colonnes du fichier
	f1B.close()
	
	for line in allLines:
		CHROM = line.split("	")[0]
		POS1 = line.split("	")[1] 
		POS2 = line.split("	")[2]
		IN_FILE = line.split("	")[3] 
		REF1 = line.split("	")[4]
		REF2 = line.split("	")[5]   
		ALT1 = line.split("	")[6] 
		ALT2 = line.split("	")[7] 
		if IN_FILE == "B":
			f1 = open(pathF1,"a") # Ouverture du fichier "No Stop"
			f1.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %(CHROM,POS1,POS2,IN_FILE,REF1,REF2,ALT1,ALT2)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
			f1.close()
		elif IN_FILE != "B":
			f1B = open(pathF1B,"a") # Ouverture du fichier "No Stop"
			f1B.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %(CHROM,POS1,POS2,IN_FILE,REF1,REF2,ALT1,ALT2)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
			f1B.close()

	################# Fichier C ####################
	
	c_sites = open(pathC,"r") 
	cLines = c_sites.readlines() 
	c_sites.close()

	f2 = open(pathF2,"w") # Création d'un nouveau fichier
	f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  | brin \n' %('chromosome','old_start','old_end','new_start','new_end','codon','Length add')) # Nomme les colonnes du fichier
	f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))
	
	f2B = open(pathF1B,"w") # Création d'un nouveau fichier
	f2B.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %('CHROM','POS1','POS2','IN_FILE','REF1','REF2','ALT1','ALT2')) # Nomme les colonnes du fichier
	f2B.close()
	
	for line in cLines:
		CHROM = line.split("	")[0]
		POS1 = line.split("	")[1] 
		POS2 = line.split("	")[2]
		IN_FILE = line.split("	")[3] 
		REF1 = line.split("	")[4]
		REF2 = line.split("	")[5]   
		ALT1 = line.split("	")[6] 
		ALT2 = line.split("	")[7] 
		if IN_FILE == "B":
			f2 = open(pathF1,"a") # Ouverture du fichier "No Stop"
			f2.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %(CHROM,POS1,POS2,IN_FILE,REF1,REF2,ALT1,ALT2)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
			f1.close()
		elif IN_FILE != "B":
			f2B = open(pathF1B,"a") # Ouverture du fichier "No Stop"
			f2B.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %(CHROM,POS1,POS2,IN_FILE,REF1,REF2,ALT1,ALT2)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
			f2B.close()

	################# Fichier J ####################
	
	j_sites = open(pathJ,"r") 
	jLines = j_sites.readlines() 
	j_sites.close()

	f3 = open(pathF3,"w") # Création d'un nouveau fichier
	f3.write('%15s | %15s | %15s | %15s | %15s | %15s  | brin \n' %('chromosome','start','end','first_frame','last_frame','brin')) # Nomme les colonnes du fichier
	f3.write('%15s | %15s | %15s | %15s | %15s | %15s  | ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))
	f3B = open(pathF1B,"w") # Création d'un nouveau fichier
	f3B.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %('CHROM','POS1','POS2','IN_FILE','REF1','REF2','ALT1','ALT2')) # Nomme les colonnes du fichier
	f3B.close()
	
	for line in jLines:
		CHROM = line.split("	")[0]
		POS1 = line.split("	")[1] 
		POS2 = line.split("	")[2]
		IN_FILE = line.split("	")[3] 
		REF1 = line.split("	")[4]
		REF2 = line.split("	")[5]   
		ALT1 = line.split("	")[6] 
		ALT2 = line.split("	")[7] 
		if IN_FILE == "B":
			f3 = open(pathF1,"a") # Ouverture du fichier "No Stop"
			f3.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %(CHROM,POS1,POS2,IN_FILE,REF1,REF2,ALT1,ALT2)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
			f3.close()
		elif IN_FILE != "B":
			f3B = open(pathF1B,"a") # Ouverture du fichier "No Stop"
			f3B.write('%15s%15s%15s%15s%15s%15s%15s%15s \n' %(CHROM,POS1,POS2,IN_FILE,REF1,REF2,ALT1,ALT2)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
			f3B.close()