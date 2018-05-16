## Written by Elise GUERET in 2018


## This WDL pipeline implements removing sequencing duplicates and base recalibration
## derived from the Best Practices for Variant Discovery in RNAseq (Mai 2015) 
## and adapted to GATK 4 for preparing Dicentrarchus labrax RNA-seq data for variant analysis.


Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

# Workflow Definition

workflow Data_PreProcessing {
  
  File picard
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  File refFasta
  File refIndex
  File refDict
  File variationSites

  scatter (sample in inputSamples) {
    call ReorderSam {
      input:
        sampleName=sample[0],
        bamFile=sample[1], 
        bamIndex=sample[2], 
        RefFasta=refFasta, 
        PICARD=picard,
        RefIndex=refIndex ,
        RefDict=refDict
    }
    call Markduplicates {
      input: 
        sampleName=sample[0],
        RefFasta=refFasta, 
        PICARD=picard, 
        Reordereds=ReorderSam.Reordered
    }
    call SortSAM {
      input: 
        sampleName=sample[0], 
        RefFasta=refFasta, 
        PICARD=picard,
        MarkDups=Markduplicates.MarkDup 
        
    }
  }
}

# Task Definition

# Mets dans le même ordre les fichiers BAM que le fichier fasta du génome de référence.
task ReorderSam {
  File PICARD
  File bamFile
  File bamIndex
  File RefFasta
  File RefIndex
  String sampleName
  File RefDict
  command {
    java -jar ${PICARD} \
      ReorderSam \
      I=${bamFile} \
      O=${sampleName}_reordered.bam \
      R=${RefFasta} \
      CREATE_INDEX=TRUE
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File Reordered = "${sampleName}_reordered.bam"
  }
}
   
#
task Markduplicates {
  File PICARD
  Array[File] Reordereds
  File RefFasta
  String sampleName
  command {
    java -jar ${PICARD} \
      MarkDuplicates \
      R=${RefFasta} \
      I=${sep= "I=" Reordereds} \
      CREATE_INDEX=true \
      REMOVE_SEQUENCING_DUPLICATES=true	\
      O=${sampleName}_marked_duplicates.bam \
      M=${sampleName}_marked_dup_metrics.txt
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File MarkDup = "${sampleName}_marked_duplicates.bam"
    File Metrics = "${sampleName}_marked_dup_metrics.txt"
  }
}

task SortSAM {
  File PICARD
  File RefFasta
  Array[File] MarkDups
  String sampleName
  command {
    java -jar ${PICARD} \
      SortSam \
      R=${RefFasta} \
      I=${sep= "I=" MarkDups} \
      O=${sampleName}_marked_duplicates_sorted.bam \
      SORT_ORDER=coordinate 
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File BamSorted = "${sampleName}_marked_duplicates_sorted.bam"
  }
}




