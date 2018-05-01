# Compseq
import sys
import argparse
from Fonction import motif,sequence


def command():
	''' Cette fonction permet de prendre en compte les différentes options du programme'''
	parser = argparse.ArgumentParser(description='Calculate the composition of unique words in sequences')
	parser.add_argument('-sequence', action="store", dest='seq',type=str, required=True, help='Sequence(s) filename and optional format')
	parser.add_argument('-word', action="store", dest='nb', type=int, required=True, help='[2] This is the size of word (n-mer) to \ncount. \nThus if you want to count codon frequencies \nfor a nucleotide sequence, you should enter \n3 here. (Integer 1 or more)')
	parser.add_argument('-frame', action="store", dest='fr', type=int, default = 0, help="[0] The normal behaviour of 'compseq' is to count the frequencies of all words that occur by moving a window of length 'word' up by one each time. This option allows you to move the window up by the length of the word each time, skipping over the intervening words. You can count only those words that occur in a single frame of the word by setting this value to a number other than zero. If you set it to 1 it will only count the words in frame 1, 2 will only count the words in frame 2 and so on. (Integer 0 or more)")
	parser.add_argument('-ignorebz', action="store", dest='bz',type=str,default = 'True' , help="[Y] The amino acid code B represents Asparagine or Aspartic acid and the code Z represents Glutamine or Glutamic acid. These are not commonly used codes and you may wish not to count words containing them, just noting them in the count of 'Other' words.")
	parser.add_argument('-reverse', action="store", dest='re', type=bool,default = False,  help="[N] Set this to be true if you also wish to also count words in the reverse complement of a nucleic sequence.")
	parser.add_argument('-calcfreq', action="store", dest='cf', type=bool,default = False, help="[N] If this is set true then the expected frequencies of words are calculated from the observed frequency of single bases or residues in the sequences. If you are reporting a word size of 1 (single bases or residues) then there is no point in using this option because the calculated expected frequency will be equal to the observed frequency. Calculating the expected frequencies like this will give an approximation of the expected frequencies that you might get by using an input file of frequencies produced by a previous run of this program. If an input file of expected word frequencies has been specified then the values from that file will be used instead of this calculation of expected frequency from the sequence, even if 'calcfreq' is set to be true.")
	parser.add_argument('-zerocount', action="store", dest='zc',type=str, default = 'True', help="[Y] You can make the output results file much smaller if you do not display the words with a zero count.")
	parser.add_argument('-sbegin1', action="store", dest='begin',type=int, default = 1, help='[1]Start of each sequence to be used')
	parser.add_argument('-send1', action="store", dest='end',type=int, default = 1, help='[1]End of each sequence to be used')
	parser.add_argument('-scircular1', action="store", dest='ci',type=bool,default = False, help='[Y]Sequence is circular')
	parser.add_argument('-stdout', action="store", dest='sd',type=bool,default = False, help='[Y]Write first file to standard output')
	parser.add_argument('-version', action="store_true", dest='vs',default = False, help='EMBOSS:6.6.0.0')
	parser.add_argument('-sreverse1', action="store", dest='sr1', type=bool,default = False,  help="Reverse (if DNA)")
	return parser

def clacfreq(seq,largeur):
	''' Permet de calculer la fréquence de chaque résidu dans une séquence et de les ajouter dans un dictionnaire'''
	seq = seq.upper()
	result = {}	
	base = 'ABCDEFGHIKLMNPQRSTUVWYZ' # Comprend les sigles pour les Acides aminés mais aussi les acides nucléiques de l'ADN
	base = base.upper()
	for AA in base :
		result[AA] = 0
		for elt in seq :
			if elt == AA:
				result[AA]+=1
		result[AA] = result[AA]/largeur
	return result # Renvoie le dictionnaire de fréquence de chaque résidu
	
	
def tableau(fichier,n,reverse,sbegin1,frame,send1,circular,calcfreq,zerocount,ignorebz,sreverse1):
	''' Permet d'afficher un tableau donnant les fréquences observées et attendues de tous les motifs possibles'''
	tableau = '#'+'\n%-7s %-15s %-15s %-15s %-15s' % ('# Word','Obs Count','Obs Frequency','Exp Frequency','Obs/Exp Frequency') +'\n#'+'\n'
	other = 0
	seq = ''
	seqRev = ''
	largeur = 0
	longueur = 0
	l = len(sequence(fichier,n,reverse,sbegin1,send1,circular,sreverse1))
	for i in sequence(fichier,n,reverse,sbegin1,send1,circular,sreverse1) : # Boucle qui permet d'ajouter toutes les séquences pour le calcul lorsque le fichier fasta en contient plusieurs
		a, Seq,SeqRev = i
		seq = seq +'\n' + Seq
		if reverse == True :
			seqRev = seqRev +'\n'+ SeqRev
		else :
			seqRev = SeqRev
		longueur = (len(Seq) -(frame))//n + longueur # Permet de connaitre le nombre d'analyse effectuées si l'option frame est activée
		largeur = largeur + len(Seq)
	if frame == 0 :
		longueur = (largeur - (l*(n-1))) # Permet d'avoir le nombre de motif analysé lorsque frame est désactivée
	if frame != 0 :
		step = n
	else :
		step = 1 # Lorsque l'option frame n'est pas activée, par défaut le pas est de 1
	for i in motif(n,seq,ignorebz):
		nbr = 0
		for j in range(frame+1,largeur,step):
			if seq[j:j+n] == i :
				nbr += 1  	
				other += 1		 # Nombre de répétition du motif trouvé dans la séquence
		for k in range(frame +1,len(seqRev)-(n),step):	
			if seqRev[k:k+n] == i :
				nbr += 1  	
				other += 1
		if reverse == True :
			fobs = nbr/(longueur*2) # Calcul de la fréquence observée
		else :
			fobs = nbr/(longueur)
		if calcfreq == True :
			if reverse == True :
				seq1 = seq + '\n' + seqRev
				P = clacfreq(seq1,largeur*2) # P est le dictionnaire de fréquence pour chaque résidu
			else :
				P = clacfreq(seq,largeur) # P est le dictionnaire de fréquence pour chaque résidu
			fth = 1
			for AA in i : # Ici AA représente chaque résidus du motif choisis
				fth = fth * P[AA] # On multiplie les fréquences de chaque résidu du motif
		else :
			if 'S' and 'R' and 'L' in seq :
				fth = 1/(21**n) # Calcul de la fréquence attendue pour une protéine
			else :
				fth = 1/(4**n) # Calcul de la fréquence attendue pour une sequence nucléique
		if fth == 0 :
			ratio = 10000000000.0000000
		else :
			ratio = fobs/fth # Ratio de la fréquence observée sur la fréquence attendue
		if zerocount == 'False':
			if nbr != 0 :
				tableau = tableau + '%-7s %-15d %-15.7f %-15.7f %-15.7f' % (i,nbr,fobs,fth,ratio)+'\n'
		else :
			tableau = tableau + '%-7s %-15d %-15.7f %-15.7f %-15.7f' % (i,nbr,fobs,fth,ratio)+'\n'
	noth =longueur- other # Nombre d'autres motifs non pris en compte ailleurs dans le tableau
	fobs = noth/longueur # calcul de la fréquence observée pour d'autres motifs non pris en compte ailleurs dans le tableau
	fth =  noth//3/longueur # Calcul de la fréquence attendue pour d'autres motifs non pris en compte ailleurs dans le tableau
	ratio = 10000000000.0000000 # ratio de la fréquence observée sur la fréquence attendue
	tableau = tableau + '\n'+'%-7s %-15d %-15.7f %-15.7f %-15.7f'%('Other',noth,fobs,fth,ratio) + '\n'
	return tableau	
	
def entete(fichier,n,reverse,sbegin1,frame,send1,circular,calcfreq,zerocount,ignorebz,sreverse1):
	''' Donne l'entete du tableau obtenu par la fonction resultat'''
	id= ''
	seq = ''
	seqRev = ''
	largeur = 0
	longueur = 0
	l = len(sequence(fichier,n,reverse,sbegin1,send1,circular,sreverse1))
	for i in sequence(fichier,n,reverse,sbegin1,send1,circular,sreverse1):
		ide, Seq,SeqRev = i
		longueur = (len(Seq) -(frame))//n + longueur # Permet de connaitre le nombre d'analyse effectuées si l'option frame est activée
		largeur = largeur + len(Seq)
		id = id + '\n#	' + ide # Permet d'avoir les id de chaque séquence dans l'entête du fichier résultat
		seq = Seq +'\n' + seq  # Permet d'additionner les séquences pour l'analyse lorsque le fichier en contient plusieurs
	
	if frame == 0 :
		longueur = (largeur - (l*(n-1))) # Permet de connaitre le nombre d'analyse effectuées
	if reverse == True :
		longueur = longueur *2  # Permet de connaitre le nombre d'analyse effectuer si le brin complémentaire est ajouté à l'analyse
	a="\n# The Expected frequencies are calculated on the (false) assumption that every \n# word has equal frequency. \n# \n# The input sequences are:	" + id + '\n \n \nWord size	' + str(n) + '\nTotal count	' + str(longueur) + '\n'
	if calcfreq == True : # Permet de changer l'entête pour indiquer que l'option calcfreq est activée
		a="\n# The Expected frequencies are calculated from the observed single\n# base or residue frequencies in these sequences \n# \n# The input sequences are: 	" + id + '\n \n \nWord size	' + str(n) + '\nTotal count	' + str(longueur) + '\n'
	if frame !=0 :# Permet de changer l'entête pour indiquer que l'option frame est activée
		a = '\n# Only words in frame '+str(frame)+' will be counted.' +a
	if ignorebz == 'False' :# Permet de changer l'entête pour indiquer que l'option ignorebz est activée
		a = "\n# The amino acid codes 'B' and 'Z' will be counted,"+ "\n# rather than treated as 'Other'." + a
	if zerocount == 'False':# Permet de changer l'entête pour indiquer que l'option zerocount est activée
		a = '\n# Words with a frequency of zero are not reported.' + a
	a = "# \n# Output from 'compseq'\n#"+ a
	return a
	
def resultat():
	'''Envoie tous les résultats dans un fichier txt'''
	parser=command()
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	Arguments= parser.parse_args() # Permet, pour la suite, de donner une variable à l'argument pris en compte par la fonction command()
	fichier = Arguments.seq
	n = Arguments.nb
	reverse = Arguments.re
	sbegin1 =Arguments.begin
	send1 = Arguments.end
	frame = Arguments.fr
	circular = Arguments.ci
	calcfreq = Arguments.cf
	zerocount = Arguments.zc
	stdout = Arguments.sd
	ignorebz = Arguments.bz
	version = Arguments.vs
	sreverse1 = Arguments.sr1
	id , seq,seqRev = sequence(fichier,n,reverse,sbegin1,send1,circular,sreverse1)[0] # On récupère seulement l'ID de la première séquence pour le nom de notre fichier
	if 'S' and 'R' and 'L' in seq : # Permet d'enlever les options qui ne sont pas nécessaires pour des séquences protéiques
		sreverse1 = False
		circular = False
		reverse = False
		if n >= 4 : # Permet de limiter la recherche a un motif de 3 résidus et d'afficher un message d'avertissement si cette valeur est supérieure ou égale à 4
			print('Warning: integer value out of range ' +str(n)+ ' more than (reset to) 3')
			n =3
	else :
		if n >= 7: # Permet de limiter la recherche a un motif de 6 résidus et d'afficher un message d'avertissement si cette valeur est supérieure ou égale à 7
			print('Warning: integer value out of range '+str(n)+' more than (reset to) 6')
			n=6
	if frame >= (n +1) : # Permet de limiter l'option frame a la longueur du motif et d'afficher un message d'avertissement si cette valeur est supérieure ou égale à la longueur du motif
		print('Warning: integer value out of range '+str(frame)+' more than (reset to) '+ str(n))
		frame = n
	if version == True : # Permet d'afficher la version du logiciel
		print('EMBOSS:version Elise et Flo')
	if stdout == False : # Permet de donner les résultats dans le fichier de sortie
		fout = id +'.composition' # Permet de donner un nom au fichier de sortie avec l'id de la séquence + .composition 
		f = open(fout,'w')
		f.write(entete(fichier,n,reverse,sbegin1,frame,send1,circular,calcfreq,zerocount,ignorebz,sreverse1)+ '\n' + str(tableau(fichier,n,reverse,sbegin1,frame,send1,circular,calcfreq,zerocount,ignorebz,sreverse1))) # Permet d'écrire dans le fichier de sortie la fonction entete puis la fonction Tableau
		f.close() 	
	else : # Permet d'afficher le résultat dans le terminal et non dans le fichier de sortie
		print(entete(fichier,n,reverse,sbegin1,frame,send1,circular,calcfreq,zerocount,ignorebz,sreverse1)+ '\n' + str(tableau(fichier,n,reverse,sbegin1,frame,send1,circular,calcfreq,zerocount,ignorebz,sreverse1)))
		
			 
if __name__ == '__main__':	
	resultat()
