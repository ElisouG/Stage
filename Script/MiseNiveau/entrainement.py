#!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package entrainement.py
# @author Elise GUERET


################# String type #####################
a = "je m'appelle Elise"
a = a.upper() # Mets tout en majuscule
a = a.lower() # Mets tous en minuscule
b = 15
c = a + str(b) # Change un chiffre en lettre
d = "et j'ai"
e = "ans"
print(a+" "+d+" "+str(b)+" "+e)
print('%s %s %d %s' % (a,d,b,e))
text = a[len(a)-5:len(a)]
print(text)
fileName = 'BRCA1.fasta'
fileName = fileName.replace('.fasta','')
print(fileName)

################## List type ###################

liste = ["et","j'habite","a","montpellier"]
liste.sort()
#print(liste)
liste.append("depuis")
print(liste)
liste2 = ['pas','tres tres','longtemps']
liste.extend(liste2)
liste[6] = 'tres'
print(liste[0:8])
