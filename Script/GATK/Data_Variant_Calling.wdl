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
		call select as selectSNPs {
			input:
				sampleName=sample[0], 
				RefFasta=refFasta, 
				GATK=gatk, 
				RefIndex=refIndex, 
				RefDict=refDict, 
				type="SNP",
				rawVCFs=HaplotypeCallerERC.rawVCF
		}
		call hardFilterSNP {
			input: 
				sampleName=sample[0], 
				RefFasta=refFasta, 
				GATK=gatk, 
				RefIndex=refIndex, 
				RefDict=refDict, 
				rawSNPs=selectSNPs.rawSubset
		}
		call select as selectIndels {
			input: 
				sampleName=sample[0], 
				RefFasta=refFasta, 
				GATK=gatk, 
				RefIndex=refIndex, 
				RefDict=refDict, 
				type="INDEL", 
				rawVCFs=HaplotypeCallerERC.rawVCF
		}
		call hardFilterIndel {
			input: 
				sampleName=sample[0], 
				RefFasta=refFasta, 
				GATK=gatk, 
				RefIndex=refIndex, 
				RefDict=refDict, 
				rawIndels=selectIndels.rawSubset
		}
		call VariantsToTable {
  			input:
  	  			sampleName=sample[0],
  	  			RefFasta=refFasta,
  	  			GATK=gatk,
  	  			rawSNPs=selectSNPs.rawSubset
  		}
  	}
	call CombineGVCFs {
		input:
			RefFasta=refFasta, 
			GATK=gatk, 
			RefIndex=refIndex, 
			RefDict=refDict,
			filteredIndels=hardFilterIndel.filteredIndel,
			filteredSNPs=hardFilterSNP.filteredSNP
	}
	call GenotypeGVCFs {
		input: 
			GATK=gatk, 
			RefFasta=refFasta, 
			RefIndex=refIndex, 
			RefDict=refDict,
			VariantCombineds=CombineGVCFs.VariantCombined
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

#This task calls GATK's tool, SelectVariants, in order to separate indel calls from SNPs in
#the raw variant vcf produced by HaplotypeCaller. The type can be set to "INDEL"
#or "SNP".
task select {
	File GATK
	File RefFasta
	File RefIndex
	File RefDict
	String sampleName
	String type
	Array[File] rawVCFs
	command {
		java -jar ${GATK} \
			SelectVariants \
			-R ${RefFasta} \
			-V ${sep="-V" rawVCFs} \
			-select-type ${type} \
			-O ${sampleName}_raw.${type}.vcf
	}
	output {
		File rawSubset = "${sampleName}_raw.${type}.vcf"
	}
}
    
task VariantsToTable {
	File GATK
	File RefFasta
  	String sampleName
	Array[File] rawSNPs

	command {
	java -jar ${GATK} \
	  VariantsToTable \
      -V ${sep="-V" rawSNPs} \
      -O ${sampleName}.snps.indels.table
	}
	output {
	File TableFilteredVCF = "${sampleName}.snps.indels.table"
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
	Array[File] rawSNPs
	command {
		java -jar ${GATK} \
			VariantFiltration \
			-V ${sep="-V" rawSNPs} \
			--filter-expression "FS > 60.0" \
			--filter-name "snp_filter" \
			-O ${sampleName}.filtered.snps.vcf
	}
	output {
		File filteredSNP = "${sampleName}.filtered.snps.vcf"
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
	Array[File] rawIndels
	command {
		java -jar ${GATK} \
			VariantFiltration \
			-V ${sep="-V" rawIndels} \
			--filter-expression "FS > 200.0" \
			--filter-name "indel_filter" \
			-O ${sampleName}.filtered.indels.vcf
	}
	output {
		File filteredIndel = "${sampleName}.filtered.indels.vcf"
	}
}

# This task calls GATK's tool, .
task CombineGVCFs {
	File GATK
	File RefFasta
	File RefIndex
	File RefDict
	Array[File] filteredIndels
	Array[File] filteredSNPs
	command {
		java -jar ${GATK} \
			CombineGVCFs \
			-V ${sep="-V" filteredIndels} \
			-V ${sep="-V" filteredSNPs} \
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