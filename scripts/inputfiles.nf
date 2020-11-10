reads = Channel.fromPath( '/Users/subinpark/fingerprinting/D0/*.fastq.gz' )

process foo {
    input:
    file sample from reads.collect()
    script:
    """
    echo your_command --reads $sample
    """
}