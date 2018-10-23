#!/home/egueret/miniconda3/bin/python3
##!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Modif-CDS_start_end.py
# @author Elise GUERET


#qsub -cwd -V -S /bin/bash -l h_rt=99:00:00 -M elise.gueret@gmail.com -m bes -N Test_sge_Script_Correction -o /home/egueret/Stage_UM_ISEM/sge_test_Script_Correction.out -e /home/egueret/Stage_UM_ISEM/sge_test_Script_Correction.err -b y "python3 /home/egueret/Stage/Script/Modif-CDS_start_end.py"

if __name__ == "__main__":

#################### Path File Cluster    ######################

	#pathAnnotation = "/home/egueret/Stage_UM_ISEM/Puce_57K/correction_annotation.gtf"
	pathAnnotation = "/mnt/c/Users/missl/Documents/correction_annotation.gtf"

################### Récupération info #################

	info = open(pathAnnotation,"r")
	lines = info.readlines() # Sauvegarder dans la variable lines les lignes du fichier précédent
	info.close() # Fermer le fichier

	listeFct = [] # Création d'une liste vide
	listeFinale = [] # Création d'une liste vide
	listeStartTrue = [] # Création d'une liste vide
	listeStartFalse = [] # Création d'une liste vide
	listeStopTrue = [] # Création d'une liste vide
	listeStopFalse = [] # Création d'une liste vide
	listeAnomalies = [] # Création d'une liste vide
	listeCorrigee = [] # Création d'une liste vide
	start = False # Initialisation du compteur
	petitGene = False
	exonAfter = False
	nbExon = 0 
	nbCDS = 0

	for line in lines:
		lineSplit = line.split('\t')
		K = lineSplit[0]
		brin = lineSplit[6]
		gene_id = lineSplit[8].split(" ")[1].replace('"',"")
		fct = lineSplit[2]
		st = lineSplit[3]
		end = lineSplit[4]
		listeFct.append([gene_id,brin,K,fct,st,end])
		if brin == '+' :
			# Rechercher les 1ers CDS, Start et exon
			if fct == 'start_codon' :
				# listeFinale.append([gene_id,brin,K,fct,st,end])
				startCodon = [st,end]
				firstCDS = firstExon =	lastCDS = lastExon = ''
				nbExon = nbCDS = 0
				start = True
			elif fct == 'CDS' and start == True :
				firstCDS = [st,end]
				# listeFinale.append([gene_id,brin,K,fct,st,end])
			elif fct == 'exon' and start == True :
				firstExon = [st,end]
				# listeFinale.append([gene_id,brin,K,fct,st,end])
				start = False
			# Rechercher les derniers CDS et exon + le stop
			elif fct == 'CDS' :
				nbCDS += 1
				lastCDS =  [st,end]
			elif fct == 'exon' and petitGene == True :
				firstExon = [st,end]
				listeFinale.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				petitGene = False
				# listeFinale.append(lastExon)
				# listeFinale.append(lastCDS)
				# listeFinale.append([gene_id,brin,K,fct,st,end])
			elif fct == 'exon' and exonAfter == True and old_Gene_id == gene_id :
				lastExon = [st,end]
				listeFinale.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				exonAfter = False
			elif fct == 'exon' :
				nbExon += 1 
				lastExon = [st,end]
			elif fct == 'stop_codon' :
				stopCodon = [st,end]
				if firstExon == '' :
					petitGene = True 
				elif nbExon != nbCDS or firstCDS == '' :
					exonAfter = True
					old_Gene_id = gene_id
				else :
					listeFinale.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])



		if brin == "-":
			# Rechercher les derniers CDS et exon + le stop
			if fct == 'stop_codon' :
				stopCodon = [st,end]
				# listeFinale.append([gene_id,brin,K,fct,st,end])
				firstCDS = 	firstExon = lastCDS = lastExon = ''
				nbExon = nbCDS = 0
				start = True
			elif fct == 'CDS' and start == True :
				lastCDS = [st,end]
				nbCDS += 1
				# listeFinale.append([gene_id,brin,K,fct,st,end])
			elif fct == 'exon' and start == True :
				lastExon = [st,end]
				nbExon += 1
				# listeFinale.append([gene_id,brin,K,fct,st,end])
				start = False
			# Rechercher les 1ers CDS, Start et exon
			elif fct == 'CDS' :
				nbCDS += 1
				firstCDS = [st,end]
			elif petitGene == True and fct == 'exon':
				lastExon = [st,end]
				listeFinale.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				petitGene = False
			elif fct == 'exon' and exonAfter == True and old_Gene_id == gene_id :
				firstExon = [st,end]
				listeFinale.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				exonAfter = False
			elif fct == 'exon' :
				nbExon += 1
				firstExon = [st,end]
			elif fct == 'start_codon' :
				startCodon = [st,end]
				if lastExon == '' :
					petitGene = True
				elif nbExon != nbCDS or lastCDS == '' :
					exonAfter = True
					old_Gene_id = gene_id
				else :
					listeFinale.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])

######### Test des starts ##############

	for elt in listeFinale :
		if '' in elt :
			listeAnomalies.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
			#print('Anomalies',elt)
		else :	
			gene_id = elt[0]
			brin = elt[1]
			startCodon = elt[2]
			firstCDS = elt[3]
			firstExon = elt[4] 
			lastCDS = elt[5] 
			lastExon = elt[6]
			stopCodon = elt[7]
			if brin == '+' :
				if startCodon[0] == firstExon[0] == firstCDS[0] :
					listeStartTrue.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
					#print('True(+)',elt)
				elif startCodon[0] != firstExon[0] or firstExon[0] != firstCDS[0] or startCodon[0] !=  firstCDS[0] : 
					listeStartFalse.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
					#print('False(+)',elt)
	
			elif brin == '-' :
				if startCodon[1] == firstExon[1] == firstCDS[1] :
					listeStartTrue.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
					#print('True(-)',elt)
				elif startCodon[1] != firstExon[1] or firstExon[1] != firstCDS[1] or startCodon[1] !=  firstCDS[1] :
					listeStartFalse.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
					#print('False(-)',elt)

################ Test des Falses ###################

	for elt in listeFinale :
		gene_id = elt[0]
		brin = elt[1]
		startCodon = elt[2]
		firstCDS = elt[3]
		firstExon = elt[4] 
		lastCDS = elt[5] 
		lastExon = elt[6]
		stopCodon = elt[7]
		print('stopCodon',stopCodon)
		print('lastExon',lastExon)
		print('lastCDS',lastCDS)
		if brin == '+' :
			if stopCodon[0] == lastExon[0] == lastCDS[0] :
				listeStopTrue.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				#print('True(+)',elt)
			elif stopCodon[0] != lastExon[0] or lastCDS[0] != lastCDS[0] or stopCodon[0] !=  lastCDS[0] :
				listeStopFalse.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				#print('False(+)',elt)
		elif brin == '-' :
			if stopCodon[1] == lastExon[1] == lastCDS[1] :
				listeStopTrue.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				#print('True(+)',elt)
			elif stopCodon[1] != lastExon[1] or lastCDS[1] != lastCDS[1] or stopCodon[1] !=  lastCDS[1] : 
				listeStopFalse.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				#print('False(+)',elt)


################### CORRECTION #################

	for elt in listeStartFalse :
		gene_id = elt[0]
		brin = elt[1]
		startCodon = elt[2]
		firstCDS = elt[3]
		firstExon = elt[4] 
		lastCDS = elt[5] 
		lastExon = elt[6]
		stopCodon = elt[7]
		if brin == '+' :
			#if startCodon[0] != firstExon[0] or firstExon[0] != firstCDS[0] or startCodon[0] !=  firstCDS[0] :
			startCodon[0] = firstExon[0] = firstCDS[0]
			listeCorrigee.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
			#print('Corrigee(+)',elt)
		elif brin == '-' :
			#if startCodon[1] != firstExon[1] or firstExon[1] != firstCDS[1] or startCodon[1] !=  firstCDS[1] :
			startCodon[1] = firstExon[1] = firstCDS[1]
			listeCorrigee.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
			#print('Corrigee(-)',elt)

################### VERIFICATION #################

	for elt in listeCorrigee :
		gene_id = elt[0]
		brin = elt[1]
		startCodon = elt[2]
		firstCDS = elt[3]
		firstExon = elt[4] 
		lastCDS = elt[5] 
		lastExon = elt[6]
		stopCodon = elt[7]
		if brin == '+' :
			if startCodon[0] == firstExon[0] == firstCDS[0] :
				listeStartTrue.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				#print('True(+)',elt)
			else :
				listeStartFalse.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				#print('False(+)',elt)
		elif brin == '-' :
			if startCodon[1] == firstExon[1] == firstCDS[1] :
				listeStartTrue.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				#print('True(-)',elt)
			else :
				listeStartFalse.append([gene_id,brin,startCodon,firstCDS,firstExon,lastCDS,lastExon,stopCodon])
				#print('False(-)',elt)

	# for elt in listeStopFa:
	# 	print(elt)



######## Règle générale des conditions ###############################

# if A and B --> A et B doivent être vraie

# if A or B --> A ou B doit être vraie

# if A an B or C --> (A et B) ou C doivent être vraie

# if A and B or A and C 
# ou 
# if A  :
#	if B or C :

##############################################################


	# 	if gene_id == gene_id_prev :
	# 		gene_id = gene_id
	# 		brin = brin
	# 		k = K 	
	# 		if brin == '+' :
	# 			if fct == 'start_codon' :
	# 				listeFinale.append()
	# 			if fct == "CDS" or fct == "exon" :
	# 				fct = lineSplit[2]
	# 				st = lineSplit[3]
	# 				end = lineSplit[4]
	# 		elif brin == '-' and fct == "CDS" or fct == "exon" and a == True:
	# 			a = False
	# 			fct = lineSplit[2]
	# 			st = lineSplit[3]
	# 			end = lineSplit[4]

	# 	K_prev = K 
	# 	brin_prev = brin
	# 	gene_id_prev = gene_id




	# 	listeFinale.append([gene_id,brin,K,fct,st,end])
	# print(listeFinale)






