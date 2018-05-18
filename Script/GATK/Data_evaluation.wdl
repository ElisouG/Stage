## Written by Elise GUERET in 2018


## This WDL pipeline implements variant calling (GVCF generation) and hardfiltering
## derived from the Best Practices for Variant Discovery in RNAseq (Mai 2015) 
## for SNP and Indel discovery in Dicentrarchus labrax RNA-seq data.

# Workflow Definition

workflow DataEvaluation {

	File picard
	File refFasta
	File refIndex
	File refDict
	File dbSNP

	call CollectVariantCallingMetrics {
		input:
			PICARD=picard,
			gVCF=GenotypeGVCFs.finalVCF,
			RefFasta=refFasta, 
			RefIndex=refIndex, 
			RefDict=refDict,
			DBSNP=dbSNP
	}
}

# Task Definition

task CollectVariantCallingMetrics {

	File PICARD
	File RefFasta
	File RefIndex
	File RefDict
	File DBSNP

	command {
		java -jar ${PICARD} \
			CollectVariantCallingMetrics \
			NSNP=${DBSNP} \
			I=${gVCF} \
			O=Variants_final_metrics
	}
	output {
		File 
	}
}
 task GenotypeConcordance {

 	File PICARD

 	command {
 		java -jar ${PICARD} \
 			GenotypeConcordance \
 			CALL_VCF=input.vcf \ 
 			CALL_SAMPLE=sample_name \ 
 			O=Variants_final_gc \ 
 			TRUTH_VCF=truth_set.vcf \ 
 			TRUTH_SAMPLE=sample_in_truth \ 
 			INTERVALS=confident.interval_list \ 
 			MISSING_SITES_HOM_REF = true
 	}
 	output {
 		File 
 	}
}