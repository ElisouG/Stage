#!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package fonctions.py
# @author Elise GUERET

import sys
import random


def seqfasta(fichier):
	''' permet d'extraire la/les séquence(s) d'un fichier fasta'''
	fichier = fichier.strip()
	fichier = open(fichier,"r")
	resultat = []
	toutesLignes = fichier.readlines()
	fichier.close()
	x= 0
	seq=""
	for uneLigne in toutesLignes[0:] : # Boucle qui permet de créer une liste de tuple contenant l'id et la séquence
		if uneLigne[0] == '>' :
			if x == 0:
				ide = uneLigne.rstrip("\n").strip(">")
		elif uneLigne != '\n':
			x+=1
			seq = seq + uneLigne.strip("\n")
	resultat.append([ide,seq])
	return resultat

def Reverse(seq) :
	''' Permet de faire le brin compélementaire de la séquence double brin'''
	seq = seq.upper()
	seqRev = seq.replace('A','t').replace('T','a').replace('C','g').replace('G','c')
	seqRev = seqRev.upper()
	seqRev1 = ''
	for i in range(1,len(seqRev)+1):
		seqRev1 =  seqRev1 + seqRev[-i] 
	return seqRev1

def sequence(fichier,n,reverse,sbegin1,send1,circular,sreverse1): 
	''' Fonction qui permet de traiter la séquence en fonction de son format 
	et de la modifier en fonction des demandes ( reverse, adn circulaire...)'''
	if 'fasta' in fichier : # Condition qui permet de traiter les fichiers fasta
		sequences = seqfasta(fichier)
	else :
		sequences = seqembl(fichier)
	Resultat = []
	sbegin1 = sbegin1 -1
	for i in sequences: # Cette boucle permet de traiter chaque tuple de la liste (id, séquence) un par un
			id,seq = i
			seq = seq.upper()
			if sreverse1 == True :
				seq = Reverse(seq)
			if circular == True :
				seq = seq + seq[0:n-1] # Permet de rajouter les résidus nécessaires (pour le motif ou le cadre de lecture) du début de la séquence à la fin de la séquence
			if send1 != 1 :
				seq = seq[0:send1] # Permet de finir la séquence au résidu indiqué par l'option send1
			seq = seq.replace(seq[:sbegin1],'',sbegin1) # Permet de commencer au résidu indiqué par l'option sbegin1
			seqRev = ''
			if reverse == True :
				seqRev = Reverse(seq)# Permet de faire la séquence complémentaire de notre séquence
			Resultat.append((id,seq,seqRev)) # Créer une liste de tuple à 3 éléments contenant l'id, la séquence et la séquence complémentaire	
	return Resultat