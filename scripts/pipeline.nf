params.reads='/Users/subinpark/fingerprinting/D0/*_{R1,R2}_001.fastq.gz'
#reads_ch=Channel.fromFilePairs(params.reads).view()
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir="results_test1"

println "reads: $params.reads"

log.info """\
           
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()