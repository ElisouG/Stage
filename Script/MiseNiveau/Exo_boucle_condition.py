##!/usr/bin/python3
#!/root/miniconda3/bin/python3
# -*- coding: utf-8 -*-
# @package Exo_boucle_condition.py
# @author Florian Charriat

################ Exercice 1, les boucles ################

# 1) Crée une liste 1 qui contient toutes les listes ci-dessous en éléments

listeA = ['AG0004','+','165004','166544','Genome1']
listeB = ['CH1857','+','58964789','58966598','Genome1']
listeC = ['GY11','-','7589574','7590685','Genome2']
listeD = ['IR009','+','1256896','1257635','Genome2']
listeE = ['FR13-myc','-','256863','257008','Genome2']
listeF = ['GN0006','+','34568','35980','Genome1']
listeG = ['GN0001','-','45896542','45897458','Genome3']
listeH = ['ZN0006','+','16564789','16565698','Genome1']

liste1 = [listeA,listeB,listeC,listeD,listeE,listeF,listeG,listeH]

# 2) Crées 5 liste : id, brin, start,end, genome
ID = []
brin = []
start = []
end = []
genome = []
for elt in liste1:
	ID.append(elt[0])
	brin.append(elt[1])
	start.append(elt[2])
	end.append(elt[3])
	genome.append(elt[4])


# 3) Utilise la première liste pour imprimer un tableau des liste suivantes (astuce voir les anciens script avec les %)
print('%20s | %20s | %20s | %20s | %20s' %('ID','brin','start','end','genome'))
print('%20s | %20s | %20s | %20s | %20s' %('-'*20,'-'*20,'-'*20,'-'*20,'-'*20))
for elt in liste1:
	print('%20s | %20s | %20s | %20s | %20s' %(elt[0],elt[1],elt[2],elt[3],elt[4]))


# 4) Utilise les 5 autres listes crée dans la question 2 pour faire la même chose que dans la questionj 3 ( astuce boucle avec 'i')
print('%20s | %20s | %20s | %20s | %20s' %('ID','brin','start','end','genome'))
print('%20s | %20s | %20s | %20s | %20s' %('-'*20,'-'*20,'-'*20,'-'*20,'-'*20))
for i in range(0,len(ID)-1):
	print('%20s | %20s | %20s | %20s | %20s' %(ID[i],brin[i],start[i],end[i],genome[i]))

# 5) Utilise les boucles pour me donner la longueur de chaque gène

length = []
for elt in liste1:
	elt[3] = int(elt[3])
	elt[2] = int(elt[2])
	elt.append(elt[3]-elt[2])
print('%20s | %20s | %20s | %20s | %20s | %20s' %('ID','brin','start','end','genome','length'))
print('%20s | %20s | %20s | %20s | %20s | %20s' %('-'*20,'-'*20,'-'*20,'-'*20,'-'*20,'-'*20))
for elt in liste1:
	print('%20s | %20s | %20s | %20s | %20s | %20s' %(elt[0],elt[1],elt[2],elt[3],elt[4],elt[5]))
# 6) Voir avec flo : trié les liste en fonction de la position start

listeSort = sorted(liste1,key=lambda start: start[2])
print('%20s | %20s | %20s | %20s | %20s | %20s' %('ID','brin','start','end','genome','length'))
print('%20s | %20s | %20s | %20s | %20s | %20s' %('-'*20,'-'*20,'-'*20,'-'*20,'-'*20,'-'*20))
for elt in listeSort:
	print('%20s | %20s | %20s | %20s | %20s | %20s' %(elt[0],elt[1],elt[2],elt[3],elt[4],elt[5]))

# 7) Compter le nombre de genome2 

nbGenome2 = genome.count('Genome2')
print(nbGenome2)

##################### Les conditions ########################

### Pour l'affichage dans cette partie tu peux te servir du code réalisé dans la question 3 de l'exercice boucle

# 1) Utilise les liste ci dessus pour afficher la liste seulement si c'est le genome1 et sinon afficher : 'Mauvaise genome'

for elt in liste1:
	if 'Genome1' in elt:
		print(elt)
	else :
		print("Mauvais Génome")

# 2) Utilise les liste ci dessus pour affichier la liste seulement si c'est le genome 1 et le brin +

for elt in liste1:
	if 'Genome1' in elt and '+' in elt:
		print (elt)


# 3) Utilise les liste ci dessus pour affichier la liste seulement si c'est  le brin + et le genome 1 ou le genome2

for elt in liste1:
	if 'Genome1' in elt or 'Genome2' in elt and '+' in elt:
		print (elt)

# 4) Affiche en premier les brins + puis en second les brin -

liste1Sorted = sorted(liste1,key=lambda brin: brin[1])
for elt in liste1Sorted:
	print(elt)

# nbCDS = 0
# nbStart =0
# lines = f.readlines()
# for line in lines:
# 	nbCDS = nbCDS + line.counts('CDS')
# 	nbStart = nbStart + line.count('start')

