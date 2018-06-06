## Written by Elise GUERET in 2018


## This WDL pipeline implements removing sequencing duplicates and base recalibration
## derived from the Best Practices for Variant Discovery in RNAseq (Mai 2015) 
## and adapted to GATK 4 for preparing Dicentrarchus labrax RNA-seq data for variant analysis.


Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

# Workflow Definition

workflow Data_PreProcessing {
  
  File gatk
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
        RefDict=refDict,
        RefIndex=refIndex, 
        PICARD=picard, 
        Reordereds=ReorderSam.Reordered
    }
    call SortSAM {
      input: 
        sampleName=sample[0], 
        RefFasta=refFasta,
        RefDict=refDict,
        RefIndex=refIndex, 
        PICARD=picard,
        MarkDups=Markduplicates.MarkDup 
        
    }
    call BaseRecalibrator {
      input: 
        sampleName=sample[0], 
        RefFasta=refFasta,
        RefDict=refDict,
        RefIndex=refIndex, 
        GATK=gatk,
        VariationSites=variationSites,
        BamSorteds=SortSAM.BamSorted
        
    }
    call ApplyBQSR {
      input: 
        sampleName=sample[0], 
        RefFasta=refFasta,
        RefDict=refDict,
        RefIndex=refIndex, 
        BaseRecals=BaseRecalibrator.BaseRecal,
        BamSorteds=SortSAM.BamSorted,
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

# This task calls Picard's tool, ReoderSam. This tool classes reads in the bam file in the same 
# order of the reference genome. 
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
   
# This task calls Picard's tool, Markduplicates. This tool marks duplicates and remove 
# sequencing duplicates from a bam file to another bam file.
task Markduplicates {
  File PICARD
  Array[File] Reordereds
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  command {
    java -Xmx50g -jar ${PICARD} \
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

# This task calls Picard's tool, SortSam. This tool permits to classify reads by coordinate.
task SortSAM {
  File PICARD
  File RefFasta
  File RefIndex
  File RefDict
  Array[File] MarkDups
  String sampleName
  command {
    java -jar ${PICARD} \
      SortSam \
      R=${RefFasta} \
      CREATE_INDEX=true \
      I=${sep= "I=" MarkDups} \
      O=${sampleName}_marked_duplicates_sorted.bam \
      SORT_ORDER=coordinate 
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File BamSorted = "${sampleName}_marked_duplicates_sorted.bam"
  }
}

# This task calls GATK's tool, BaseRecalibrator. This tool create a table which contains
# quality value to be change because sequencing duplicates were removed.
task BaseRecalibrator {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  Array[File] BamSorteds
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

# This task calls GATK's tool, ApplyBQSR. This tool applies the recalibration from the 
# table created by BaseRecalibrator.
task ApplyBQSR {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  Array[File] BamSorteds
  Array[File] BaseRecals
  command {
    java -jar ${GATK} \
      ApplyBQSR \
      -R /home/egueret/Stage_UM_ISEM/Donnees_CRECHE/ref/labrax.fasta \
      -I ${sep= "-I" BamSorteds} \
      --bqsr-recal-file ${sep= "--bqsr-recal-file" BaseRecals} \
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


