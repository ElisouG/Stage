Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

workflow ReassignOneMappingQualityFilter {
  
  File picard
  File gatk
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  File refFasta
  File refIndex
  File refDict

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
    call PrintReads {
      input: 
        sampleName=sample[0],
        Reordereds=ReorderSam.Reordered,
        RefFasta=refFasta, 
        GATK=gatk,
        RefIndex=refIndex ,
        RefDict=refDict
    }    
  } 
}

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
task PrintReads {
  File GATK
  File RefFasta
  File RefIndex
  Array[File] Reordereds
  String sampleName
  File RefDict
  command {
    java -jar ${GATK} \
      -T PrintReads \
      -I ${sep="-I" Reordereds} \
      -rf ReassignOneMappingQuality \
      -R ${RefFasta} \
      -o ${sampleName}_reassigned.bam \
      -RMQF 50 \
      -RMQT 60
  }
  # runtime { sge_queue: "cemeb20.q" }
  output {
    File Reassigned = "${sampleName}_reassigned.bam"
  }
}