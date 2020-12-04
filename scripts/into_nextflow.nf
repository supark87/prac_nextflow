params.reads='/Users/subinpark/fingerprinting/D0/*_{R1,R2}_001.fastq.gz'
//params.reads='/Users/subinpark/NeST_CDC/NeST/input2/*_{R1,R2}_001.fastq.gz'

params.references='/Users/subinpark/fingerprinting/reference_genes/cpmp_NC_004325.fasta'
Channel
    .fromFilePairs(params.reads)
    .into {reads_ch; reads_ch2}

//Channel
//   .fromPath(params.cpmp)
    
params.outdir="result_tset_k127"

//type, variable, channel for input and output

process fastqc {

    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_out"
    
    input:
    tuple val(sample_id), file(reads_file) from reads_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}
process bbduk{

    tag "$sample_id"
    publishDir "${params.outdir}/bbduk_out", mode:'copy'

    input:
    tuple val(sample_id), file(reads_file) from reads_ch2

    output:
    tuple val(sample_id), file("${sample_id}_clean_{1,2}.fastq") into skesa_input
    file "${sample_id}.stats.txt"

    script:
    """
    bbduk.sh -Xmx1g in=${reads_file[0]} in2=${reads_file[1]} out=${sample_id}_clean_1.fastq out2=${sample_id}_clean_2.fastq qtrim=rl trimq=30 minlength=100\
    stats=${sample_id}.stats.txt
    
    """

}


process skesa{

    tag "$sampe_id"
    publishDir "${params.outdir}/skesa_out", mode:'copy'

    input:
    tuple val(sample_id),file(reads_file) from skesa_input

    output:
    tuple val(sample_id), file("${sample_id}_skesa.fasta") into skesa_output,skesa_output_for_prod
     
    script:
    """
    skesa --reads ${reads_file[0]},${reads_file[1]} --use_paired_ends --kmer 31  --contigs_out ${sample_id}_skesa.fasta
    """
}

// process run_prodigal{
//     tag "$sample_id"
//     publishDir "${params.outdir}/prod_fasta",mode:'copy'

//     input:
//     set val(sample_id), file(assembly) from skesa_output_for_prod

//     output:
//     set val(sample_id), file("${assembly}_prod.fasta") into prod_fasta

//     when:
//     assembly.size() > 0

//     script:
//     """
//     prodigal -i ${assembly} -p meta -d ${assembly}_prod.fasta
//     """
// }

// process makeblastdb{

//     tag "$sample_id"
//     publishDir "${params.outdir}/blast_db",mode:'copy'

//     input:
//     set val(sample_id),file(assembly) from prod_fasta

//     output:
//     val "${sample_id}_db" into blastdb_name
//     path "${sample_id}_dir" into blastdb_dir

//     when:
//     assembly.size() > 0
    
//     script:
//     """
//     mkdir ${sample_id}_dir
//     makeblastdb -in ${assembly} -dbtype nucl -out ${sample_id}_dir/${sample_id}_db
//     """


// }

// process blastrun2{
//     tag "$sample_id"
//     publishDir "${params.outdir}/blast_out_files2",mode:'copy'
    
//     input:
//     path cpmp from params.references
//     val sample_id_db from blastdb_name 
//     path dbdir from blastdb_dir

//     output:
//     file ("${sample_id_db}.blast")

//     script:

//     """
//     blastn -query $cpmp -db $dbdir/$sample_id_db  -outfmt "6 qseqid sseqid pident qlen qstart qend sstart send sseq" -out ${sample_id_db}.blast
//     """
 

// }

// // process filteroutput{

// // }

// // process alginment{

// // }

// // process assign_alleles{

// // }

//value : sampleid, file: paired files, and channel
//when you run, you can put "-process.echo"and it will show the process
//-with-docker

//docker run $HOME:$HOME --workdir $PWD -u $(id -u):$(id -g) <command>

/*docker image made
docker tag <image name> supark87/<imagename:version>
docker push supark87/<imagename:version>
*/
//nextflow run ./scripts/clean.nf -with-docker supark87/my-image:v2