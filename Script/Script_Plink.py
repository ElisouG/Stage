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

# Ram to use
# -l mf=30G

# Inform scheduler how many Ram are used
# -l h_vmem=30G

# Choose to run on a specific queue
# -q cemeb20.q

# Job hardware time limit
#$ -l h_rt=99:00:00

# Job name
#$ -N "run_PLINK"

# Standard output (you can use the option '-j y' which add stderr with sdtout)
#$ -o /home/egueret/Stage_UM_ISEM/Donnees_CRECHE/sge_PLINK.out

# Erreur output (don't use with the option '-j y')
#$ -e /home/egueret/Stage_UM_ISEM/Donnees_CRECHE/sge_PLINK.err

# Analyse des fréquences alléliques
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --freq --out Freq_all

# Analyse des SNPs manquants
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --missing --out Missing

# Suppression des génotype ayant un callrate inférieur à 95%
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --make-bed -- mind 0.05 --out Highgeno

# Analyse de stratifictation de la population
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --cluster --mc2 --ppc 0.05 --out Strat

# Ancestry
#/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --mh --within pop.phe --adjust --out Ancestry

# Visualisation de la substructure
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --cluster --matrix --out Substructure

# TDT : test disease trait
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --parentdt2 --out TDT2
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --parentdt1 --out TDT1
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --tdt --out TDT

# POO : Parent of origin
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --tdt --poo --out POO

# QFAM : Association trait quantitatif (family based)
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --qfam-total --mperm 100000 --out QFAM_tot
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --qfam --mperm 100000 --out QFAM
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --qfam-parents --mperm 100000 --out QFAM_parents
/home/egueret/tools/plink-1.07-x86_64/plink --file home/egueret/Stage_UM_ISEM/vcftools/Variants_final --qfam-between --mperm 100000 --out QFAM_between
