#!/bin/bash


# Toutes
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Association_All --pheno PHENO.txt --phenoname PHENO --select --out Bayesien
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Association_All --pheno PHENO.txt --phenoname PHENO --phenops --out Frequentist
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile All_Freq --pheno All_Pheno.txt --phenoname Lignee --phenops --out Frequentist_Lignee
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile All_Freq --pheno All_Pheno.txt --phenoname sexe --phenops --out Frequentist_sexe
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile All_Freq --pheno All_Pheno.txt --phenoname EAI_class --phenops --out Frequentist_eai
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile All_Freq --pheno All_Pheno.txt --phenoname Comportement --phenops --out Frequentist_comportement

# C
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno PHENO_ligneeC_qual.txt --phenoname PHENO --select --out C_Bayesien
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno PHENO_ligneeC_qual.txt --phenoname PHENO --phenops --out C_Frequentist
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno C_All_Pheno.txt --phenoname Lignee --phenops --out C_Frequentist_Lignee
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno C_All_Pheno.txt --phenoname sexe --phenops --out C_Frequentist_sexe
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno C_All_Pheno.txt --phenoname EAI_class --phenops --out C_Frequentist_eai
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno C_All_Pheno.txt --phenoname Comportement --phenops --out C_Frequentist_comportement

# J
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno PHENO_ligneeJ_qual.txt --phenoname PHENO --select --out J_Bayesien
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno PHENO_ligneeJ_qual.txt --phenoname PHENO --phenops --out J_Frequentist
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno J_All_Pheno.txt --phenoname Lignee --phenops --out J_Frequentist_Lignee
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno J_All_Pheno.txt --phenoname sexe --phenops --out J_Frequentist_sexe
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno J_All_Pheno.txt --phenoname EAI_class --phenops --out J_Frequentist_eai
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno J_All_Pheno.txt --phenoname Comportement --phenops --out J_Frequentist_comportement

# Cpos
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno PHENO_ligneeC_qual.txt --phenoname PHENO --select --out C_Bayesien
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno PHENO_ligneeC_qual.txt --phenoname PHENO --phenops --out C_Frequentist
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Cpos_eaimoy_Association_Linear --pheno C_All_Pheno.txt --phenoname Lignee --phenops --out Cpos_Frequentist_Lignee
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Cpos_eaimoy_Association_Linear --pheno C_All_Pheno.txt --phenoname sexe --phenops --out Cpos_Frequentist_sexe
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Cpos_eaimoy_Association_Linear --pheno C_All_Pheno.txt --phenoname EAI_class --phenops --out Cpos_Frequentist_eai
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Cpos_eaimoy_Association_Linear --pheno C_All_Pheno.txt --phenoname Comportement --phenops --out Cpos_Frequentist_comportement

# Jpos
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno PHENO_ligneeJ_qual.txt --phenoname PHENO --select --out J_Bayesien
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno PHENO_ligneeJ_qual.txt --phenoname PHENO --phenops --out J_Frequentist
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Jpos_eaimoy_Association_Linear --pheno J_All_Pheno.txt --phenoname Lignee --phenops --out Jpos_Frequentist_Lignee
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Jpos_eaimoy_Association_Linear --pheno J_All_Pheno.txt --phenoname sexe --phenops --out Jpos_Frequentist_sexe
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Jpos_eaimoy_Association_Linear --pheno J_All_Pheno.txt --phenoname EAI_class --phenops --out Jpos_Frequentist_eai
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Jpos_eaimoy_Association_Linear --pheno J_All_Pheno.txt --phenoname Comportement --phenops --out Jpos_Frequentist_comportement

# Cneg
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno PHENO_ligneeC_qual.txt --phenoname PHENO --select --out C_Bayesien
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile C_Association_All --pheno PHENO_ligneeC_qual.txt --phenoname PHENO --phenops --out C_Frequentist
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Cneg_eaimoy_Association_Linear --pheno C_All_Pheno.txt --phenoname Lignee --phenops --out Cneg_Frequentist_Lignee
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Cneg_eaimoy_Association_Linear --pheno C_All_Pheno.txt --phenoname sexe --phenops --out Cneg_Frequentist_sexe
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Cneg_eaimoy_Association_Linear --pheno C_All_Pheno.txt --phenoname EAI_class --phenops --out Cneg_Frequentist_eai
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Cneg_eaimoy_Association_Linear --pheno C_All_Pheno.txt --phenoname Comportement --phenops --out Cneg_Frequentist_comportement

# Jneg
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno PHENO_ligneeJ_qual.txt --phenoname PHENO --select --out J_Bayesien
#/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile J_Association_All --pheno PHENO_ligneeJ_qual.txt --phenoname PHENO --phenops --out J_Frequentist
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Jneg_eaimoy_Association_Linear --pheno J_All_Pheno.txt --phenoname Lignee --phenops --out Jneg_Frequentist_Lignee
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Jneg_eaimoy_Association_Linear --pheno J_All_Pheno.txt --phenoname sexe --phenops --out Jneg_Frequentist_sexe
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Jneg_eaimoy_Association_Linear --pheno J_All_Pheno.txt --phenoname EAI_class --phenops --out Jneg_Frequentist_eai
/home/elise.gueret/Téléchargements/trinculo_v0.96/bin/trinculo multinom --bfile Jneg_eaimoy_Association_Linear --pheno J_All_Pheno.txt --phenoname Comportement --phenops --out Jneg_Frequentist_comportement
