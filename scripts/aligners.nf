params.aligner = 'kallisto'

process foo {
  script:
  if( params.aligner == 'kallisto' )
    """
    kallisto --reads /some/data.fastq
    """
  else if( params.aligner == 'salmon' )
    """
    salmon --reads /some/data.fastq
    """
  else
    throw new IllegalArgumentException("Unknown aligner $params.aligner")
}