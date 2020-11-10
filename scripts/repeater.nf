sequences = Channel.fromPath("$baseDir/data/prots/*tfa")
methods=['regular','expresso','psicoffee']

process alignSequences {
    input:
    path seq from sequences
    each mode from methods
    // each : repeater

    """
    echo t_coffee -in $seq -mode $mode 
}