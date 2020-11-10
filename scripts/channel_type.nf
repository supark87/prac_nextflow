//que cahannel cannot be used more than one

//ch=Channel.from(1,2,3)
//println(ch) // groovy funtion
//ch.view()
//ch.view()

//  vs. Values channels can be used for many times : good for using reference files

//ch=Channel.value('Hello')
//ch.view()
//ch.view()

//from

//ch=Channel.from(1,3,5,7)
//ch.view{"value:$it"}

// of : range of valeues

//Channel
//    .of(1..23,'X','Y')
//    .view()

//list

//fromPath

//Channel.fromPath('/Users/subinpark/fingerprinting/D0/*_{R1,R2}_001.fastq.gz',hidden:true) 
 
 
 //   .view()


//reads_ch=Channel.fromFilePairs('/Users/subinpark/fingerprinting/D0/*_{R1,R2}_001.fastq.gz',
//,flat: true, checkIfExists:true)
//                .view()

// flat : no list any more




