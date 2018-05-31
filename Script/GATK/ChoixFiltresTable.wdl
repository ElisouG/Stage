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
		call VariantsToTable {
  			input:
  	  			sampleName=sample[0],
  	  			RefFasta=refFasta,
  	  			GATK=gatk,
  	  			rawSNPs=selectSNPs.rawSubset
  		}
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
      -F CHROM -F POS -F QUAL -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum -F InbreedingCoeff  \
      -O ${sampleName}.snps.indels.table     
	}
	output {
	File TableFilteredVCF = "${sampleName}.snps.indels.table"
	}
}

