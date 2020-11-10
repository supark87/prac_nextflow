/*
 * pipeline input parameters
 */
params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

println """\
         R N A S E Q - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


/*
 * define the `index` process that create a binary index
 * given the transcriptome file
 */
process index {
    //cpus 2 

    input:
    path transcriptome from params.transcriptome

    output:
    path 'index' into index_ch

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

//task.cpus: modify via configs file
//if you run all 
// nextflow run script.nf --reads 'data/ggal/*_{1,2}.fq'

read_pairs_ch=Channel.fromFilePairs(params.reads,checkIfExists:true)
/*or use operater
Channel.fromFilePairs(params.reads)
       .set{read_pairs_ch}

read_paris_ch.view
/-resume : when you add reads, it will skipped already done