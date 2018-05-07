Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

workflow RemoveDuplicates {
  
  File gatk
  File picard
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  File refFasta
  File refIndex
  File refDict

  scatter (sample in inputSamples) {
    call Markduplicates {
      input: 
        sampleName=sample[0],
        bamFile=sample[1], 
        bamIndex=sample[2], 
        RefFasta=refFasta, 
        PICARD=picard, 
    }    
  } 
}
    
#
task Markduplicates {
  File PICARD
  File bamFile
  File bamIndex
  File RefFasta
  String sampleName
  command {
    java -jar ${PICARD} \
      MarkDuplicates \
      R=${RefFasta} \
      I=${bamFile} \
      CREATE_INDEX=true \
      REMOVE_SEQUENCING_DUPLICATES=true	\
      O=${sampleName}_marked_duplicates.bam \
      M=${sampleName}_marked_dup_metrics.txt
  }
  runtime { sge_queue: "cemeb20.q" }
  output {
    File MarkDup = "${sampleName}_marked_duplicates.bam"
    File Metrics = "${sampleName}_marked_dup_metrics.txt"
  }
}