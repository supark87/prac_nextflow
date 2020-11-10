//singularity pull docker <>

And use 

nextflow run -with-singularity <docker image>

//Conda

nextflor run <script> -with-conda <conda env>

//You can define container for each process

process A{
    container <name>

}

process B{
    container <name>
}

command -profile(calling from config file)