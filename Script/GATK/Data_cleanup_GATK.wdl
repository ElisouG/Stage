## Written by Elise GUERET in 2018


## This WDL pipeline implements base recalibration
## derived from the Best Practices for Variant Discovery in RNAseq (Mai 2015) 
## and adapted to GATK 4 for preparing Dicentrarchus labrax RNA-seq data for variant analysis.


Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

# Workflow Definition

workflow DataCleanupGATK {
  
  File gatk
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  File refFasta
  File refIndex
  File refDict
  File variationSites

  scatter (sample in inputSamples) {
    call BaseRecalibrator {
      input: 
        sampleName=sample[0], 
        RefFasta=refFasta, 
        GATK=gatk,
        VariationSites=variationSites,
        BamSorteds=sample[1],
        BamIndex=sample[2]
        
    }
    call ApplyBQSR {
      input: 
        sampleName=sample[0], 
        RefFasta=refFasta, 
        BaseRecals=BaseRecalibrator.BaseRecal,
        BamSorteds=sample[1],
        BamIndex=sample[2],
        GATK=gatk
        
    }
    #call AnalyseCovariate {
    	#input:
    	  #sampleName=sample[0],
        #BaseRecals=BaseRecalibrator.BaseRecal,
    	  #GATK=gatk
    #}
  }
}

# Task Definition

task BaseRecalibrator {
  File GATK
  File RefFasta
  String sampleName
  Array[File] BamSorteds
  Array[File] BamIndex
  File VariationSites
  command {
    java -jar ${GATK} \
      BaseRecalibrator \
      -I ${sep="-I" BamSorteds} \
      -R /home/egueret/Stage_UM_ISEM/Donnees_CRECHE/ref/labrax.fasta \
      -OBI true \
      --known-sites /home/egueret/Stage_UM_ISEM/Donnees_CRECHE/inputs/Final_list_57907_SNPs.recode.vcf \
      -O ${sampleName}_marked_duplicates_sorted_recal_data.table
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File BaseRecal = "${sampleName}_marked_duplicates_sorted_recal_data.table"
  }
}

task ApplyBQSR {
  File GATK
  File RefFasta
  String sampleName
  Array[File] BamSorteds
  Array[File] BamIndex
  Array[File] BaseRecals
  command {
    java -jar ${GATK} \
      ApplyBQSR \
      -R /home/egueret/Stage_UM_ISEM/Donnees_CRECHE/ref/labrax.fasta \
      -I ${sep="-I" BamSorteds} \
      --bqsr-recal-file ${sep="--bqsr-recal-file" BaseRecals} \
      -O ${sampleName}_marked_duplicates_sorted_recalibrated.bam
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File BamRecal = "${sampleName}_marked_duplicates_sorted_recalibrated.bam"
  }
}

#task AnalyseCovariate {
 # File GATK
  #Array[File] BaseRecals
  #String sampleName
  #command {
  #  java -jar ${GATK} \
  #    AnalyzeCovariates \
  #    -bqsr ${sep="-bqsr" BaseRecals} \
  #    -plots ${sampleName}_recalibration.pdf \
  #    -csv ${sampleName}_recalibration.csv
  #}
  # runtime { sge_queue: "cemeb20.q" }
  #output {
  #  File Plot = "${sampleName}_recalibration.pdf"
  #  File Csv = "${sampleName}_recalibration.csv"
  #}
#}


