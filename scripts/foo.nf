

process foo {
    //container=foo
    label 'small'
    echo true
    """
    echo Using $task.cpus cpus and $task.memory memory
    """  
}


process bar{
    //scontainer=bar
    label 'big'
    echo true
    """
    echo using $task.cpus cpus and $task.memory memory
    """
   
}
 