## Written by Elise GUERET in 2018


# Workflow Definition

workflow SplitVCF_singleSample {
  
  File gatk
  File refFasta
  File refIndex
  File refDict
  String sampleName
  File vcfFile

  call SelectVariants {
    input:
      RefFasta=refFasta, 
      GATK=gatk,
      RefIndex=refIndex,
      RefDict=refDict,
      vcfFile=vcfFile,
      sampleName=sampleName
  }
}

# Task Definition

# This task calls Picard's tool, ReoderSam. This tool classes reads in the bam file in the same 
# order of the reference genome. 
task SelectVariants {
  File GATK
  File vcfFile
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  command {
    java -jar ${GATK} \
      SelectVariants \
      -I ${vcfFile} \
      -O ${sampleName}_tobefiltered.vcf \
      -sn ${sampleName} \
      -R ${RefFasta} 
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File VCF = "${sampleName}_tobefiltered.vcf"
  }
}
   



