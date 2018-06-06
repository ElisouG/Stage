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
		call VariantsToTable {
  			input:
  	  			sampleName=sample[0],
  	  			RefFasta=refFasta,
  	  			GATK=gatk,
  	  			rawVCFs=HaplotypeCallerERC.rawVCF
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
			-stand-call-conf 10.0 \
			--dont-use-soft-clipped-bases true \
			-I ${bamFile} \
			-O ${sampleName}_rawLikelihoods.g.vcf 
	}
	output {
		File rawVCF = "${sampleName}_rawLikelihoods.g.vcf"
	}
}

task VariantsToTable {
	File GATK
	File RefFasta
  	String sampleName
	Array[File] rawVCFs
	command {
	java -jar ${GATK} \
	  VariantsToTable \
      -V ${sep="-V" rawVCFs} \
      -F CHROM -F POS -F REF -F ALT -F TYPE -F QUAL -F DP -F MLEAC -F MLEAF -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum -F InbreedingCoeff -F GT -F GQ -F MIN_DP -F PL -F SB -F RAW_MQ -F AD -F ExcessHet -F BaseQRankSum -F ClippingRankSum \
      -O ${sampleName}.snps.indels.table     
	}
	output {
	File TableFilteredVCF = "${sampleName}.snps.indels.table"
	}
}

