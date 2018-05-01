#!/usr/bin/python3.4
# -*- coding: utf-8 -*-
# @package Groupe.py
# @author Florian Charriat

from module_elise import form
################ Exercice 1, les listes ################

# 1) Crée deux listes avec au moins 3 éléments. Chaque liste doit avoir le même nombre d'element

liste1 = ["et", "il", "pleut", "toujours"]
liste2 = ["je", "travaille", "au", "bureau"]

# 2) Crée une troisième liste qui contiendras deux élements, la liste1 et la liste2 puis une quatrième liste qui contiendra tous les elements des liste 1 et 2.

liste3 = [liste1] + [liste2]
#print(liste3)

liste4 = liste1 + liste2
#print(liste4)

# 3) De la liste 3, récupère dans une variable le premier element de chaque liste de la liste.

elt1 = liste3[0][0]
elt2 = liste3[1][0]

# 4) La dernière liste est une liste contenant toutes les lettres de la variable 2

Variable2 = 'Coucou mon coeur, travaille bien'

liste5= list(Variable2)
liste5.remove(' ')
liste5.remove(' ')
liste5.remove(' ')
liste5.remove(' ')
liste5.remove(',')

# print(liste5)

################ Exercice 2 Les strings ###########################


line1 = 'id : AG0004 | frame = 3, position : 150100-180021 | chromosome : 5 , genome : M.oryzae'
line2 = 'id : BRCA2 | frame = 0, position : 10-5890 | chromosome : 12 , genome : O.sativa'
line3 = 'id : AG0004 | frame = 3, position : 150100000-1800000021 | chromosome : 2 , genome : homme'

# 1) Recupère seulement l'id, la position (sous la forme chromosome:position) et le genome.

idLine1 = line1.split("|")[0].replace("id : ", "").strip()
posLine1 =  line1.split(":")[3].split(",")[0].strip() + ":" +line1.split(":")[2].split("|")[0].strip()
genome1 = line1.split(":")[-1].strip()

idLine2 = line2.split("|")[0].replace("id : ", "").strip()
posLine2 =  line2.split(":")[3].split(",")[0].strip() + ":" +line2.split(":")[2].split("|")[0].strip()
genome2 = line2.split(":")[-1].strip()

idLine3 = line3.split("|")[0].replace("id : ", "").strip()
posLine3 =  line3.split(":")[3].split(",")[0].strip() + ":" +line3.split(":")[2].split("|")[0].strip()
genome3 = line3.split(":")[-1].strip()

# 2) Affiche sous forme de "tableau" les données recupéré avec un titre en entete de chaque colonne

aff = "%8s | %30s | %15s\n%8s | %30s | %15s\n%8s | %30s | %15s\n%8s | %30s | %15s" % ("ID","Position","Genome",idLine1, posLine1, genome1,idLine2, posLine2, genome2,idLine3, posLine3, genome3)
print(form(aff,type = ["bold",'underline'],col="blue"))
