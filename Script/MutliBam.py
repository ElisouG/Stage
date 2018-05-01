import os
from module_elise import createDir, form

pathBam = 'path/to/directory/'
pathOutput = 'path/to/directory/output'

createDir(pathOutput)
f = open('%sscript.sh'%pathOutput,'w')
f.close()
for file in os.listdir(pathBam):
	if file.endswith('.bam'):
		pathFile = '%s%s'%(pathBam,file)
		f = open('%sscript.sh'%pathOutput,'a')
		f.write('qsub -q cemeb.q ... "samttols -f %s -o %s.bai"\n'%(pathFile,pathOutput+file.replace('.bam','')))
		f.close()


print(form('Veuillez taper la commance suivante pour lancer les jobs ;',type='bold'))
print(form('\tbash %sscript.sh'%pathOutput),'green','bold')

#Commande a taper pour scp (quand tu es dans le fichier ou se trouve le dossier a envoyé)
# cmd : scp -r dossier egueret@162.38.180.18:/homedir/egueret/....