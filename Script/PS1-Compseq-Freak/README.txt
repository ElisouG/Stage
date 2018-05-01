Compseq et Freak sont des modules du logiciel EMBOSS.
Ces modules existent au format C et sont désormais en partie codés au format python3.

Le module Compseq permet d'obtenir un fichier contenant les différentes fréquences des différents motifs d'une séquence d'ADN ou d'une séquence protéique.
Ce module prend en compte différentes fonctions qui sont :
	- sequence : la séquence d'ADN ou protéique au format FASTA ou EMBL
	- word : est un motif de 1 à 6 résidus à rechercher dans la séquence entrée. Ce motif ne peut dépasser 6 pour une séquence d'ADN et 3 pour une séquence protéique.
	- scircular1 : circularise la séquence entrée
	- send1 : détermine le résidu où s'arrête la recherche du motif
	- sbegin1 : détermine le résidu où commence la recherche du motif
	- reverse : permet d'obtenir la séquence complémentaire d'ADN et l'ajoute à la séquence entrée
	- sreverse1 : permet d'obtenir la séquence complémentaire et de n'effectuer la recherche du motif que sur cette séquence complémentaire
	- frame : permet de décaler le pas de lecture du motif le long de la séquence avec un pas indiqué par la valeur du frame
	- calcfreq : permet de calculer les véritables fréquences d'observations de chacuns des résidus sur la séquence entrée
	- ignorebz : permet de comptabiliser dans une séquence protéique les acides aminés B et Z qui sont normalement ignorés par défaut
	- zerocount : permet d'éliminer du fichier de sortie les lignes pour lesquelles le ou les motifs ne sont pas présents dans la séquence entrée
	- stdout : permet de sortir les résultats directement dans la console et non dans un fichier de sortie

Le module Freak permet d'obtenir dans un fichier de sortie la fréquence de chaque résidu recherché dans une séquence d'ADN ou une séquence protéique.
Ce module prend en compte différentes fonctions qui sont les suivantes :
	- sequence : la séquence d'ADN ou protéique au format FASTA ou EMBL
	- window : c'est la fenêtre de lecture pour effectuer la recherche des résidus un à un
	- scircular1 : circularise la séquence entrée
	- send1 : détermine le résidu où s'arrête la recherche du motif
	- sbegin1 : détertermine le résidu où commence la recherche du motif
	- sreverse1 : permet d'obtenir la séquence complémentaire et de n'effectuer la recherche du motif que sur cette séquence reversée
	- letters : les différents résidus à analyser
	- stdout : permet de sortir les résultats directement dans la console et non dans un fichier de sortie

Limites:

Nos programmes peuvent bien entendu être optimisés pour un meilleur fonctionnement.
Notre fonction motif fonctionne rapidement avec des motifs jusqu'à 6 résidus pour les acides nucléiques.
Mais cette fonction devient lente voire très lente lorsque l'on cherche des motifs de 3 ou 4 résidus d'acide aminé c'est pour cela qu'on a limité les motifs à 3 acides aminés.
Notre option scircular1 fonctionne différement du programme EMBOSS.
L'extraction des id suffit mais n'est pas totalement identique à celle donnée par EMBOSS.
Notre fonction help est appelée par "-h". 
Notre fonction command() utilisant argparse gère mal les type booléens. En effet, il renvoie la valeur True si l'option est entrée dans le terminal quel que soit la valeur qui la suit. Donc les options de type booléen de notre programme se limite à True. Les options précédentes ne peuvent pas prendre la valeur False lorsqu'elles sont rentrées dans le programme. C'est pourquoi pour les options -ignorebz et -zerocount nous avons utilisés le type 'str'. On aurait pu faire de même pour tous les types booléens.

