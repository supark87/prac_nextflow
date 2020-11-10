Channel.from(1,2,3)
        .map{it ->[it*it]}
        .view{num,sqr -> "Square of : $num is $sqr}

        //or


Channel.from(1,2,3,4)
       .map{it -> it*it}
       .view()

// modify channel element with map and feed into channel

words_ch=Channel
                .from{'hello','world'}
                .map{ word -> [word,word.size()]}
                //now it is tuple

process greeting{
    echo true

    input tuple val(word), val(word_size) from words_ch

    script:
    """
    echo Word is $word and Size is $word_size
    """
}



///now if you have file object

files_ch=Channel.fromPath('data/*.fq)
                .map{file -> [file.name,file]}
process_greeting {
    echo true
    
    input:
    tuple val(fname), file(reads) from files_ch

    script:
    """
    echo Alignment of sample $fname with $reads
    """
    
}