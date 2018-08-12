def warningParse(pathFile):
	"""
	"""
	info = open(pathFile,"r") # Lire le fichier warning obtenu par snpEff
	lines = info.readlines() # Sauvegarder dans la variable lines les lignes du fichier précédent
	info.close() # Fermer le fichier
	a = 0 # Initialisation du compteur
	b = 0 # Initialisation du compteur
	pb = [] # Création d'une liste vide
	listeInfo = [] # Création d'une seconde liste vide

	for line in lines:
		# Récupérer le type de problème
		if line[0:8] == "WARNING:" or line[0:6] == "ERROR:" : # Permet de récupérer le début de la ligne
			# Triage des problèmes en fonction de leur origine
			a = 1 # Mettre le compteur a à 1
			if "non-last" in line.split("'")[-1] : # Récupérer les multiple stop
				pb.append("Multiple Stop") # Ajouter ce pb à la liste des pb
			elif "STOP" in line.split("'")[-1] : # Récupérer les no stop
				pb.append("No Stop") # Ajouter ce pb à la liste des pb
			elif "START" in line.split("'")[-1] : # Récupérer les no start
				pb.append("No Start") # Ajouter ce pb à la liste des pb
			elif "has length" in line.split("'")[-1] : # Récupérer les incomplets
				pb.append("Incomplet") # Ajouter ce pb à la liste des pb
		# Récupérer leur cadre de lecture
		elif "rank:" in line and b!=1:
			frame1 = line.split(',')[2].split(':')[1] # Récupérer le 1er cadre de lecture
			b=1 # Mettre le compteur b à 1
		elif b == 1 and "rank:" in line :
			frame2 = line.split(',')[2].split(':')[1] # Récupérer le dernier cadre de lecture
		# Récupérer le reste des informations
		elif a == 1 and len(pb) == 1: # Vérifier si le compteur est à 1
			a = 0 # Ré-initialisation du compteur
			lineSplit = line.split(",") # Séparation des caractères de la ligne par ","
			if "-" in lineSplit[1] : # Récupération de l'information brin anti-sens
				brin = "anti-sens"
			elif "+" in lineSplit[1] : # Récupération de l'information brin sens
				brin = "sens"
			lineSplit = lineSplit[0].split(":") # Séparation des caractères de la ligne par ":"
			chromosome = lineSplit[0] # Récupération du chromosome
			lineSplit = lineSplit[1].split("-") # Séparation des caractères de la ligne par "-"
			start = int(lineSplit[0]) # Récupération du début du gène et mise en integer (chiffre)
			end = int(lineSplit[1]) # Récupération de la fin du gène et mise en integer (chiffre)
		elif line[0:5] == '\t\tCDS' and len(pb) == 1:
			b = 0 # Ré-initialisation du compteur
			for elt in pb :
				listeInfo.append([chromosome,start,end,brin,elt,frame1,frame2]) # Création d'une liste comportant des listes des différentes info sur les problèmes rencontrés
			pb = [] # Ré-initialisation de la liste
		elif a == 1 and len(pb) != 1:
			pb = [] # Ré-initialisation de la liste
	return listeInfo

def correctionSens(listeCorrection,dico_GenomeR,dico_genome,pathF1,pathF2,pathF3,pathF4):
	"""
	"""
	f1 = open(pathF1,"w") # Création d'un nouveau fichier
	f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  | brin \n' %('chromosome','old_start','old_end','new_start','new_end','codon','Length add')) # Nomme les colonnes du fichier
	f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))
	#f1.write('%s\t%s\t%s\t%s\t%s\t%s')
	f2 = open(pathF2,"w") # Création d'un nouveau fichier
	f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  | brin \n' %('chromosome','old_start','old_end','new_start','new_end','codon','Length add')) # Nomme les colonnes du fichier
	f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))
	f3 = open(pathF3,"w") # Création d'un nouveau fichier
	f3.write('%15s | %15s | %15s | %15s | %15s | %15s  | brin \n' %('chromosome','start','end','first_frame','last_frame','brin')) # Nomme les colonnes du fichier
	f3.write('%15s | %15s | %15s | %15s | %15s | %15s  | ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15))
	f4 = open(pathF4,"w") # Création d'un nouveau fichier
	f4.write('%15s | %15s | %15s | %15s | %15s | %15s  | brin \n' %('chromosome','start','end','first_frame','last_frame','brin')) # Nomme les colonnes du fichier
	f4.write('%15s | %15s | %15s | %15s | %15s | %15s  | ---------- \n' %('-'*15,'-'*15,'-'*15,'-'*15,'-'*15,'-'*15,))

	codonStop = ['TAA','TAG','TGA'] # Liste des différents codons stop existants
	codonStart = 'ATG' # Nomme le codon start
	codonStart2 = 'AAG' # Nomme le codon start alternatif
	listeNoStop = []
	listeNoStart = []
	for elt in listeCorrection :

		chromosome, startGTF, endGTF, brin, pb ,frame1, frame2 = elt
		start = startGTF # Récupération de la position du début du gène
		start1 = start - 99 # Ajout de 100bp en amont du gène
		end = endGTF # Récupération de la position de la fin du gène
		end1 = end + 999 # Ajout de 1000bp en aval du gène
	################################################ Sens #######################################

		if brin == 'sens':
			sequence = dico_genome[chromosome].seq # Si le brin est sens alors il prend le fichier "sens"
			if pb == 'No Stop' :
				# if frame2 == '1' :
				# 	end = end - 2
				# 	end1 = end1 - 2
				# elif frame2 == '2' :
				# 	end = end - 1 
				# 	end1 = end1 -1
				for i in range(end+1,end1-2,3): # Recherche des codons stop en aval du gène
					if sequence[i:i+3] in codonStop :
						Longueur = (i+3) - (int(endGTF)+1)
						f2 = open(pathF2,"a") # Ouverture du fichier "No Stop"
						f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  +\n' %(chromosome,endGTF-1,endGTF+1,i+1,i+3,sequence[i:i+3],str(Longueur))) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
						f2.close() # Fermer le fichier
						listeNoStop.append([chromosome,str(int(endGTF)-1),str(int(endGTF)+1),str(i+1),str(i+3),'+'])
						break # S'arrête dès qu'il trouve le premier codon stop

			elif pb == 'No Start':
				st = 'NA'
				en = 'NA'
				# if frame1 == '1':
				# 	start = start +1
				# 	start1 = start1 +1
				# elif frame1 == '2':
				# 	start = start + 2
				# 	start1 = start1 + 2

				for i in range(start1,start-3,3): # Recherche d0es codons start en amont du gène
					if sequence[i:i+3] == codonStart :
						st = i
						en = i+2
				if st != 'NA':
					Longueur = (int(startGTF)+1) - int((st+1))
					f1 = open(pathF1,"a") # Ouverture du fichier "No Start"
					f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  +\n' %(chromosome,startGTF+1,startGTF+3,st+1,en+1,sequence[st:en+1],str(Longueur))) # Ecrire dans ce fichier les codons start trouvés ainsi que leur position
					f1.close() # Fermer le fichier 	
					listeNoStart.append([chromosome,str(int(startGTF)+1),str(int(startGTF)+3),str(st+1),str(en+1),'+'])
		
			elif pb == 'Multiple Stop':
				f3 = open(pathF3,"a") # Ouverture du fichier "No Stop"
				f3.write('%15s | %15s | %15s | %15s | %15s | %15s  |  +\n' %(chromosome,startGTF,endGTF,frame1,frame2,brin)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
				f3.close()

			elif pb == 'Incomplet':
				f4 = open(pathF4,"a") # Ouverture du fichier "No Stop"
				f4.write('%15s | %15s | %15s | %15s | %15s | %15s  |  +\n' %(chromosome,startGTF,endGTF,frame1,frame2,brin)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
				f4.close()


	############################# Anti sens ##############################################
			elif brin == 'anti-sens':
				sequence = dico_GenomeR[chromosome].seq # Si le brin est anti-sens alors il prend le fichier "anti-sens"
				length = len(sequence)
				if pb == 'No Stop' :
					# if frame1 == '1' :
					# 	start = start + 2
					# elif frame1 == '2' :
					# 	start = start +1
					for i in range(length-1-(start-1),length-1-(start-1002),3): # Recherche des codons stop en aval du gène
						if sequence[i:i+3] in codonStop :
							Longueur = (int(startGTF)+1) - (length-(i+1))
							f2 = open(pathF2,"a") # Ouverture du fichier "No Stop"
							f2.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  -\n' %(chromosome,startGTF +1,startGTF+3,str(length-(i+3)),str(length-(i+1)),sequence[i:i+3],str(Longueur))) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
							f2.close() # Fermer le fichier
							listeNoStop.append([chromosome,str(int(startGTF)+1),str(int(startGTF)+3),str(length-(i+3)),str(length-(i+1)),'-'])
							break # S'arrête dès qu'il trouve le premier codon stop

				elif pb == 'No Start' :
					st = 'NA'
					en = 'NA'
					# if frame2 == '1':
					# 	end = end -1
					# if frame2 == '2':
					# 	end = end -2
					for i in range(length-1-(end+99),length-1-(end +1) ,3): # Recherche des codons start en amont du gène
						if sequence[i:i+3] == codonStart :
							st = i
							en = i+2
					if st != 'NA':
						Longueur = (length-(st)) - (int(endGTF)+1)
						f1 = open(pathF1,"a") # Ouverture du fichier "No Start"
						f1.write('%15s | %15s | %15s | %15s | %15s | %15s  | %15s  |  -\n' %(chromosome,endGTF-1,endGTF+1,length-(en),length-(st),sequence[st:en+1],str(Longueur))) # Ecrire dans ce fichier les codons start trouvés ainsi que leur position
						f1.close() # Fermer le fichier 
						listeNoStart.append([chromosome,str(int(endGTF)-1),str(int(endGTF)+1),str(length-(en)),str(length-(st)),'-'])

				elif pb == 'Multiple Stop':
					f3 = open(pathF3,"a") # Ouverture du fichier "No Stop"
					f3.write('%15s | %15s | %15s | %15s | %15s | %15s |  +\n' %(chromosome,startGTF,endGTF,frame1,frame2,brin)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
					f3.close()

				elif pb == 'Incomplet':
					f4 = open(pathF4,"a") # Ouverture du fichier "No Stop"
					f4.write('%15s | %15s | %15s | %15s | %15s | %15s  | +\n' %(chromosome,startGTF,endGTF,frame1,frame2,brin)) # Ecrire dans ce fichier les codons stop trouvés ainsi que leur position
					f4.close()

	return listeNoStart,listeNoStop

def newAnnotation(output,pathAnnotation,listeNoStop,listeNoStart) :
	"""
	"""
	Annotation = open(pathAnnotation)
	linesAnnotation = Annotation.readlines()
	Annotation.close()
	# print('ListeAnnotation Done')
	# print(strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	newAnnotation = open(output, "w")
	for line in linesAnnotation:
		for elt in listeNoStop:
			if elt[0] in line and '\t'+elt[1]+'\t' in line and 'stop_codon' in line: 
				line = line.replace(elt[1],elt[3])
				line = line.replace(elt[2],elt[4])
			if  elt[5] == '+'  and  elt[0] in line:
				if 'exon' in line and '\t'+elt[2]+'\t' in line  :
					line = line.replace(elt[2],elt[4])
				elif 'CDS' in line and '\t'+str(int(elt[2])-3)+'\t' in line :
					line = line.replace(str(int(elt[2])-3),str(int(elt[4])-3))
			if  elt[5] == '-' and  elt[0] in line : 
				if 'exon' in line and '\t'+elt[1]+'\t' in line :
					line = line.replace(elt[1],elt[3])
				elif 'CDS' in line and '\t'+str(int(elt[1])+3) +'\t' in line:
					line = line.replace(str(int(elt[1])+3),str(int(elt[3])+3))
		for elt in listeNoStart:
			if elt[0] in line and '\t'+elt[1]+'\t' in line and 'start_codon' in line: 
				line = line.replace(elt[1],elt[3])
				line = line.replace(elt[2],elt[4])
			if  elt[5] == '+'  and  elt[0] in line:
				if 'exon' in line and '\t'+elt[1]+'\t' in line  :
					line = line.replace(elt[1],elt[3])
				elif 'CDS' in line and '\t'+ elt[1] +'\t' in line :
					line = line.replace(elt[1],elt[3])
			if  elt[5] == '-' and  elt[0] in line: 
				if 'exon' in line and '\t'+elt[2]+'\t' in line :
					line = line.replace(elt[2],elt[4])
				elif 'CDS' in line and '\t'+ elt[2] +'\t' in line :
					line = line.replace(elt[2],elt[4])

		newAnnotation.write(line)
	newAnnotation.close()

def transcriptomeParse(pathTranscriptome,pathTranscrits) :
	"""
	"""
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
	return listeTranscrit

def TestStopCorrection(pathStopTested,listeNoStop,listeTranscrit):
	"""
	"""
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
	return listeNoStopTested

def TestStartCorrection(pathStartTested,listeNoStart,listeTranscrit):
	"""
	"""
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
	return listeNoStartTested

def recupPosCDS(pathAnntotationNew,pathCDS) :
	"""
	"""
	NewGTF = open(pathAnntotationNew)
	linesNewGTF= NewGTF.readlines()
	NewGTF.close()

	listeCDS = []
	CDSFinaux = []

	CDS = open(pathCDS, "w")
	CDS.write("%s | %s | %s | %s | %s:%s\n" % ('Chromosome','geneID','brin','frame','CDS_start','CDS_end'))
	geneID = "none"
	for line in linesNewGTF:
		if '\tCDS\t' in line :
			if geneID == "none":
				K1 = line.split('\t')[0]
				brin = line.split('\t')[6]
				CDS_start = line.split('\t')[3]
				CDS_end = line.split('\t')[4]
				frame = line.split('\t')[7]
				geneID = line.split('"')[1].split('|')[0]
				CDS.write("%s | %s | %s | %s | %s:%s" % (K1,geneID,brin,frame,CDS_start,CDS_end)) 
				listeCDS = [K1,geneID,brin,frame,'%s:%s'% (CDS_start,CDS_end)]
			elif line.split('"')[1].split('|')[0] != geneID :
				geneID = line.split('"')[1].split('|')[0]
				K1 = line.split('\t')[0]
				brin = line.split('\t')[6]
				CDS_start = line.split('\t')[3]
				CDS_end = line.split('\t')[4]
				frame = line.split('\t')[7]
				CDS.write("\n%s | %s | %s | %s | %s:%s" % (K1,geneID,brin,frame,CDS_start,CDS_end)) 
				listeCDS = [K1,geneID,brin,frame,'%s:%s'% (CDS_start,CDS_end)] 
				CDSFinaux.append(listeCDS)
			elif line.split('"')[1].split('|')[0] == geneID :
				CDS_start = line.split('\t')[3]
				CDS_end = line.split('\t')[4]
				CDS.write(" | %s:%s " % (CDS_start,CDS_end))
				listeCDS.append('%s:%s'% (CDS_start,CDS_end))  
				geneID = line.split('"')[1].split('|')[0]
	CDS.close()

	return CDSFinaux

def recupSeqCDS(pathSequenceCDS,CDSFinaux,Genome,GenomeR):
	"""
	"""
	SequenceCDS = open(pathSequenceCDS, "w")
	SequenceCDS.write("%s | %s | %s | %s\n" % ('Chromosome','geneID','brin','Sequence'))

	CDSComplete = []
	for elt in CDSFinaux:
		seqFinale = ""
		K1 = elt[0]
		geneID = elt[1]
		brin = elt[2]
		frame = elt[3]
		positionCDS=elt[4:len(CDSFinaux)]
		for elt in positionCDS: # Est-ce qu'il ne faudrait pas que j'échange cette ligne avec celle du dessous?
			if brin == '+':
				sequence = Genome[K1].seq 
				# if frame == '1':
				# 	CDS_start = int((elt.split(':')[0]))+1
				# 	CDS_end = int(elt.split(':')[1])
				# 	seqCDS = sequence[CDS_start:CDS_end] 
				# 	seqFinale = seqFinale+seqCDS
				# 	SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
				# elif frame == '2':
				# 	CDS_start = int((elt.split(':')[0]))+2
				# 	CDS_end = int(elt.split(':')[1])
				# 	seqCDS = sequence[CDS_start:CDS_end]
				# 	seqFinale = seqFinale+seqCDS
				# 	SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
				# elif frame == '0':
				CDS_start = int(elt.split(':')[0]) -1 # Correction entre GTF et python
				CDS_end = int(elt.split(':')[1]) -1
				seqCDS = sequence[CDS_start:CDS_end]
				seqFinale = seqFinale+seqCDS
				SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
			
			elif brin == '-':
				sequence = GenomeR[K1].seq	
				# if frame == '1':
				# 	CDS_start = int((elt.split(':')[0]))+1
				# 	CDS_end = int(elt.split(':')[1])
				# 	seqCDS = sequence[CDS_start:CDS_end]
				# 	seqFinale = seqFinale+seqCDS
				# 	SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
				# elif frame == '2':
				# 	CDS_start = int((elt.split(':')[0]))+2
				# 	CDS_end = int(elt.split(':')[1])
				# 	seqCDS = sequence[CDS_start:CDS_end]
				# 	seqFinale = seqFinale+seqCDS
				# 	SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
				# elif frame == '0':
				CDS_start = int(elt.split(':')[0]) -1 # Correction entre GTF et python
				CDS_end = int(elt.split(':')[1]) -1 
				seqCDS = sequence[CDS_start:CDS_end]
				seqFinale = seqFinale+seqCDS
				SequenceCDS.write("%s | %s | %s | %s\n" % (K1,geneID,brin,seqCDS))
				seqFinale = seqFinale+seqCDS
			CDSComplete.append([K1,geneID,brin,seqFinale])	#Idente le de facon a qu'il soit dans la boucle et au meêm niveau que elit brin ='-';
	SequenceCDS.close()
	return CDSComplete
