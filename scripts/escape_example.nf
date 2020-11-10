params.data = 'le monde'

process baz {
  shell:
  '''
  X='Bonjour'
  echo $X !{params.data}
  '''
}

// ''' : real path """ with $: nf variable
//shell: you can mix and match
