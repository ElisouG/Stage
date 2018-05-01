#!/usr/bin/python3
#-*- coding: utf-8 -*-
# @package Correction_Annotation.py
# @author Elise GUERET
	
from itertools import groupby

def word_frequencies(content, blacklist):
    """
    Count the number of words in a content, excluding blacklisted terms.
    Return a generator of tuples (count, word) sorted by descending frequency.

    Example::

        >>> song = 'Ob la di ob la da "rla di da" da "da"'
        >>> for count, word in word_frequencies(song, ['di']):
        ...     print "%s %s" % (count, word)
        ...
        4 da
        2 la
        2 ob
        1 rla
    """
    sorted_words = sorted(word \
                        for word in content.lower().replace('"', '').split() \
                            if word not in blacklist)
    return ((len(list(group)), word) for word, group in groupby(sorted_words))




if __name__ == "__main__":

	pathWarning = "/mnt/c/Users/missl/Documents/Stage_UM-ISEM/Puce_57K/COMBINED_ANNOTATION_FUNCTION-2014.gtf"

	info = open(pathWarning,"r") # Lire le fichier warning obtenu par snpEff
	lines = info.readlines() # Sauvegarder dans la variable lines les ligne du fichier précédent
	info.close() # Fermer le fichier
	word_frequencies(lines, [''])