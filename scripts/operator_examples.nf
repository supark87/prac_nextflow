///Mix operator 

c1 = Channel.from( 1,2,3 )
c2 = Channel.from( 'a','b' )
c3 = Channel.from( 'z' )

c1 .mix(c2,c3).view()

//flatten 
foo = [[1,2,3],[9]]
bar = [4, 5, 6]

Channel
    .from(foo, bar)
    .flatten()
    .view()

///collect(opposite of flatten)

Channel 
    .from(1,2,3,4)
    .collect()
    .view()


///grouptuple

Channel
     .from( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
     .groupTuple()
     .view()

Channel
     .from( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
     .groupTuple(by:0)
     .view()

Channel
     .from( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
     .groupTuple(by:1)
     .view()



///join
left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right= Channel.from(['Z', 6], ['Y', 5], ['X', 4])
left.join(right).view()

//if you want 'P' there as well

left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right= Channel.from(['Z', 6], ['Y', 5], ['X', 4])
left.join(right,remainder:true).view()


//Branch(split channel)

Channel
    .from(1,2,3,40,50)
    .branch {
        small: it < 10
        large: it > 10
    }
    .set { result }

 result.small.view { "$it is small" }
 result.large.view { "$it is large" }