## Written by Elise GUERET in 2018


## This WDL pipeline implements variant calling (GVCF generation) and hardfiltering
## derived from the Best Practices for Variant Discovery in RNAseq (Mai 2015) 
## for SNP and Indel discovery in Dicentrarchus labrax RNA-seq data.

Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

# Workflow Definition

workflow DataVariantCalling {

	File gatk
	File inputSamplesFile
	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
	File refFasta
	File refIndex
	File refDict

	scatter (sample in inputSamples) {
		call HaplotypeCallerERC {
			input: 
				GATK=gatk, 
				RefFasta=refFasta, 
				RefIndex=refIndex, 
				RefDict=refDict, 
				sampleName=sample[0],
				bamFile=sample[1], 
				bamIndex=sample[2]
		}
	}
	call CombineGVCFs {
		input:
			RefFasta=refFasta, 
			GATK=gatk, 
			RefIndex=refIndex, 
			RefDict=refDict,
			rawVCFs=HaplotypeCallerERC.rawVCF
	}
	call GenotypeGVCFs {
		input: 
			GATK=gatk, 
			RefFasta=refFasta, 
			RefIndex=refIndex, 
			RefDict=refDict,
			VariantCombineds=CombineGVCFs.VariantCombined
	}
	call select as selectSNPs {
		input:
			RefFasta=refFasta, 
			GATK=gatk, 
			RefIndex=refIndex, 
			RefDict=refDict, 
			type="SNP",
			rawVCFs=CombineGVCFs.rawVCF
	}
	call hardFilterSNP {
		input: 
			RefFasta=refFasta, 
			GATK=gatk, 
			RefIndex=refIndex, 
			RefDict=refDict, 
			rawSNPs=selectSNPs.rawSubset
	}
	call select as selectIndels {
		input:  
			RefFasta=refFasta, 
			GATK=gatk, 
			RefIndex=refIndex, 
			RefDict=refDict, 
			type="INDEL", 
			rawVCFs=CombineGVCFs.rawVCF
	}
	call hardFilterIndel {
		input: 
			RefFasta=refFasta, 
			GATK=gatk, 
			RefIndex=refIndex, 
			RefDict=refDict, 
			rawIndels=selectIndels.rawSubset
	}
	call VariantsToTable {
		input:
	  		RefFasta=refFasta,
	  		GATK=gatk,
	  		filteredSNP=hardFilterSNP.filteredSNP
	}
	
}

# Task Definition

#This task calls GATK's tool, HaplotypeCaller in ERC mode. This tool takes a pre-processed 
#bam file and discovers variant sites. These raw variant calls are then written to a vcf file.
task HaplotypeCallerERC {

	File GATK
	File RefFasta
	File RefIndex
	File RefDict
	String sampleName
	File bamFile
	File bamIndex
	command {
		java -jar ${GATK} \
			HaplotypeCaller \
			-ERC GVCF \
			-R ${RefFasta} \
			-I ${bamFile} \
			-O ${sampleName}_rawLikelihoods.g.vcf 
	}
	output {
		File rawVCF = "${sampleName}_rawLikelihoods.g.vcf"
	}
}

# This task calls GATK's tool, .
task CombineGVCFs {
	File GATK
	File RefFasta
	File RefIndex
	File RefDict
	Array[File] rawVCFs
	command {
		java -jar ${GATK} \
			CombineGVCFs \
			-V ${sep="-V" rawVCFs} \
			-R ${RefFasta}\
			-O Variants_Combined.vcf
	}
	output {
		File VariantCombined = "Variants_Combined.vcf"
	}
}

# This task calls GATK's tool, GenotypeGVCFs. This tool performs joint genotyping on samples pre-called 
# by HaplotypeCaller in ERC mode from e GenomicsDB workspace created by GenomicsDBImport.
task GenotypeGVCFs {

	File GATK
	File RefFasta
	File RefIndex
	File RefDict
	File VariantCombineds

	command {
		java -Xmx4g -jar ${GATK} \
			GenotypeGVCFs \
			-R ${RefFasta} \
			-V ${VariantCombineds} \
			-O Variants_final.vcf
	}
	output {
		File finalVCF = "Variants_final.vcf"
	}
}

#This task calls GATK's tool, SelectVariants, in order to separate indel calls from SNPs in
#the raw variant vcf produced by HaplotypeCaller. The type can be set to "INDEL"
#or "SNP".
task select {
	File GATK
	File RefFasta
	File RefIndex
	File RefDict
	String type
	File finalVCF
	command {
		java -jar ${GATK} \
			SelectVariants \
			-R ${RefFasta} \
			-V ${finalVCF} \
			-select-type ${type} \
			-O raw.${type}.vcf
	}
	output {
		File rawSubset = "raw.${type}.vcf"
	}
}
    
#This task calls GATK's tool, VariantFiltration. It applies certain recommended filtering 
#thresholds to the SNP-only vcf. VariantFiltration filters out any variant that is "TRUE" 
#for any part of the filterExpression (i.e. if a variant has a QD of 1.3, it would be 
#filtered out). The variant calls remain in the file, but they are tagged as not passing.
#GATK tools downstream in the pipeline will ignore filtered calls by default
task hardFilterSNP {
	File GATK
	File RefFasta
	File RefIndex
	File RefDict
	String sampleName
	File rawSNPs
	command {
		java -jar ${GATK} \
			VariantFiltration \
			-V ${rawSNPs} \
			--filter-expression "FS > 60.0" \
			--filter-name "snp_filter" \
			-O filtered.snps.vcf
	}
	output {
		File filteredSNP = "filtered.snps.vcf"
	}
}

#As above, this task calls GATK's tool, VariantFiltration. However, this one applied filters
#meant for indels only.
# Choix des filtres: https://software.broadinstitute.org/gatk/documentation/article.php?id=1255
task hardFilterIndel {
	File GATK
	File RefFasta
	File RefIndex
	File RefDict
	String sampleName
	File rawIndels
	command {
		java -jar ${GATK} \
			VariantFiltration \
			-V ${rawIndels} \
			--filter-expression "FS > 200.0" \
			--filter-name "indel_filter" \
			-O iltered.indels.vcf
	}
	output {
		File filteredIndel = "filtered.indels.vcf"
	}
}

task VariantsToTable {
	File GATK
	File RefFasta
  	String sampleName
	File filteredSNP

	command {
	java -jar ${GATK} \
	  VariantsToTable \
      -V ${filteredSNP} \
      -F CHROM -F POS -F REF -F ALT -F TYPE -F QUAL -F DP -F MLEAC -F MLEAF -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum -F InbreedingCoeff -F GT -F GQ -F MIN_DP -F PL -F SB -F RAW_MQ -F AD -F ExcessHet -F BaseQRankSum -F ClippingRankSum \
      -O snps.indels.table
	}
	output {
	File TableFilteredVCF = "snps.table"
	}
}