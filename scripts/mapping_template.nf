params.transcriptome='$baseDir/data/transcriptome'
params.reads='$baseDir/data/*_{1,2}.fq'

reads_ch=Channel.fromFilePairs(paramns.reads)

process foo{
    input:
    path transcriptome from paramns.transcriptome
    tuple val(sample_id),path(reads) from reads_ch

    script:
    """
    echo align $sample_id against $transcriptome
    """

}