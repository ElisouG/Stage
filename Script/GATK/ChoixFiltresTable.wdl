## Written by Elise GUERET in 2018


## This WDL pipeline implements variant calling (GVCF generation) and hardfiltering
## derived from the Best Practices for Variant Discovery in RNAseq (Mai 2015) 
## for SNP and Indel discovery in Dicentrarchus labrax RNA-seq data.

Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

# Workflow Definition

workflow ChoixFiltresTable {

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
			VariantCombined=CombineGVCFs.VariantCombined
	}
	call VariantsToTable {
  		input:
  	  		RefFasta=refFasta,
  	  		GATK=gatk,
  	  		finalVCF=GenotypeGVCFs.finalVCF
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
			-stand-call-conf 10.0 \
			-bamout  ${sampleName}_rawLikelihoods.bam \
			--bam-writer-type CALLED_HAPLOTYPES \
			--output-mode EMIT_ALL_SITES \
			--dont-use-soft-clipped-bases true \
			-I ${bamFile} \
			-O ${sampleName}_rawLikelihoods.g.vcf 
	}
	output {
		File rawVCF = "${sampleName}_rawLikelihoods.g.vcf"
		File BamVCF =  "${sampleName}_rawLikelihoods.bam"
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
	File VariantCombined

	command {
		java -Xmx4g -jar ${GATK} \
			GenotypeGVCFs \
			-R ${RefFasta} \
			-V ${VariantCombined} \
			-O Variants_final.vcf
	}
	output {
		File finalVCF = "Variants_final.vcf"
	}
}

task VariantsToTable {
	File GATK
	File RefFasta
  	File finalVCF
	command {
	java -jar ${GATK} \
	  VariantsToTable \
      -V ${finalVCF} \
      -F CHROM -F POS -F REF -F ALT -F TYPE -F QUAL -F DP -F MLEAC -F MLEAF -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum -F InbreedingCoeff -F GT -F GQ -F MIN_DP -F PL -F SB -F RAW_MQ -F AD -F ExcessHet -F BaseQRankSum -F ClippingRankSum \
      -O snps_indels.table     
	}
	output {
	File TableFilteredVCF = "snps_indels.table"
	}
}

