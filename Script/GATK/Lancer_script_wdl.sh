#!/bin/bash

# Define shell to use for this job
#$ -S /bin/bash

# Export all my environment variables into runtime context
#$ -V

#$ -terse

# User to inform
#$ -M elise.gueret@gmail.com

# Get a mail when the job begins, ends, or is suspended.
#$ -m bes

# Using current working directory
#$ -cwd

# Choose to run on a specific queue
#$ -q cemeb20.q

# Job hardware time limit
#$ -l h_rt=99:00:00

# Job name
#$ -N "run_gatk_"

# Standard output (you can use the option '-j y' which add stderr with sdtout)
#$ -o /home/egueret/Stage_UM_ISEM/Données_CRECHE/sge_lancer.out

# Erreur output (don't use with the option '-j y')
#$ -e /home/egueret/Stage_UM_ISEM/Tuto-GATK/Erreur/sge_lancer.err

#Modules à charger pour que le script fonctionne
module load java1.8
module load java1.8_jre


#Lancer le pipeline
java -Dconfigfile=/home/egueret/Stage_UM_ISEM/Données_CRECHE/cromwell.conf -jar /home/egueret/tools/cromwell.jar run ReassignOneMappingQualityFilter.wdl -i ReassignOneMappingQualityFilter_inputs2.json