## Written by Elise GUERET in 2018


## This WDL pipeline implements removing sequencing duplicates and base recalibration
## derived from the Best Practices for Variant Discovery in RNAseq (Mai 2015) 
## and adapted to GATK 4 for preparing Dicentrarchus labrax RNA-seq data for variant analysis.


Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

# Workflow Definition

workflow Data_PreProcessing {
  
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
        BamSorteds=SortSAM.BamSorted
        
    }
    call ApplyBQSR {
      input: 
        sampleName=sample[0], 
        RefFasta=refFasta, 
        BaseRecals=BaseRecalibrator.BaseRecal,
        BamSorteds=SortSAM.BamSorted,
        GATK=gatk
        
    }
    call AnalyseCovariate {
    	input:
    	  sampleName=sample[0],
          BaseRecals=BaseRecalibrator.BaseRecal,
    	  GATK=gatk
    }
  }
}

# Task Definition

task BaseRecalibrator {
  File GATK
  File RefFasta
  String sampleName
  Array[File] BamSorteds
  File VariationSites
  command {
    java -jar ${GATK} \
      BaseRecalibrator \
      -I ${sep= "-I" BamSorteds} \
      -R ${RefFasta} \
      --known-sites ${VariationSites} \
      -O ${sampleName}_marked_duplicates_sorted_recal_data.table
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File BaseRecal = "${sampleName}_recal_data.table"
  }
}

task ApplyBQSR {
  File GATK
  File RefFasta
  String sampleName
  Array[File] BamSorteds
  Array[File] BaseRecals
  command {
    java -jar ${GATK} \
      ApplyBQSR \
      -R ${RefFasta} \
      -I ${sep= "-I" BamSorteds} \
      --bqsr-recal-file ${sep= "--bqsr-recal-file" BaseRecals} \
      -O ${sampleName}_marked_duplicates_sorted_recalibrated.bam
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File BamRecal = "${sampleName}_marked_duplicates_sorted_recalibrated.bam"
  }
}

task AnalyseCovariate {
  File GATK
  Array[File] BaseRecals
  String sampleName
  command {
    java -jar ${GATK} \
      AnalyzeCovariates \
      -bqsr ${sep= "--bqsr-recal-file" BaseRecals} \
      -plots ${sampleName}_recalibration.pdf  
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File Plot = "${sampleName}_recalibration.pdf"
  }
}


