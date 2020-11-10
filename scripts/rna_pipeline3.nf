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
*/

process quant {
    tag "$pair_id"
     
    input:
    file index from index_ch
    set pair_id, file(reads) from read_pairs_ch
 
    output:
    file(pair_id) into quant_ch
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

process fastqc {
    tag "FASTQC on $sample_id"

    input:
    set sample_id, file(reads) from read_pairs2_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}  
  
  
process multiqc {
    publishDir params.outdir, mode:'copy'
       
    input:
    path '*' from quant_ch.mix(fastqc_ch).collect() 

    output:
    'multiqc_report.html'
     
    script:
    """

    multiqc . 
    """
}

workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}

//(base) ➜  prac_nextflow git:(main) ✗ nextflow run nextflow-io/rnaseq-nf
//(base) ➜  prac_nextflow git:(main) ✗ nextflow run nextflow-io/rnaseq-nf -r v1.2

mail {
  from = 'info@nextflow.io'
  smtp.host = 'email-smtp.eu-west-1.amazonaws.com'
  smtp.port = 587
  smtp.user = "xxxxx"
  smtp.password = "yyyyy"
  smtp.auth = true
  smtp.starttls.enable = true
  smtp.starttls.required = true
}

//mkdir bin
//put script to bin and you can customize it

