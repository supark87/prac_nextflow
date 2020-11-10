params.reads='/Users/subinpark/fingerprinting/D0/*_{R1,R2}_001.fastq.gz'
reads_ch=Channel.fromFilePairs(params.reads).view()

//type, variable, channel for input and output

process fastq {
    
    tag "$sample_id"
    publishDir "results/fastqc", mode: 'copy'    
   
    input:
    tuple val(sample_id), file(reads_file) from reads_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    echo mkdir fastqc_${sample_id}_logs
    echo fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}

//value : sampleid, file: paired files, and channel
//when you run, you can put "-process.echo"and it will show the process