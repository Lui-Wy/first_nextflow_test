nextflow.enable.dsl = 2 

params.temp = "${projectDir}/downloads" 
params.out = "${projectDir}/output"
params.with_fastqc = false 
params.with_stats = false


process prefetch {
    //storeDir params.temp
    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"

    output:
        path "SRR1777174/SRR1777174.sra" 

    """
    prefetch SRR1777174
    """
}

process split {
    //storeDir params.temp
    publishDir params.out, mode: "copy", overwrite: true

    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"
    input:
        path inputvariable
    output:
        path "${inputvariable.getSimpleName()}.fastq"

    """
    fastq-dump --split-3 "${inputvariable}" 
    """
//split-3 does already generate an output - so I don´t need a redirection like > "test.txt"

}


process ngsutils {
    container "https://depot.galaxyproject.org/singularity/ngsutils%3A0.5.9--py27heb79e2c_4"
    //storeDir params.temp // not such a good idea? if I change things, but I run the process before - the output is the unchanged stuff in the storeDir not the new stuff!
    publishDir params.out, mode: "copy", overwrite: true
   input:
        path inputvariable
    output:
        path "${inputvariable}_stats.txt"

    """
    fastqutils stats "${inputvariable}" > "${inputvariable}_stats.txt"
    """

}

process fastqc {
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    //storeDir params.temp
    publishDir params.out, mode: "copy", overwrite: true
   input:
        path inputvariable
    output:
        path "${inputvariable.getSimpleName()}*" 

    """
    fastqc -o . -f fastq "${inputvariable}" 
    """

}



//    usage of fastqc:
// fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] [-c contaminant file] seqfile1 .. seqfileN



workflow {
    //c_prefetch = prefetch()
    //c_split = split(c_prefetch)
    //c_split.view()
    //ngsutils(c_split)
    //fastqc(c_split)


    if (params.with_fastqc == false && params.with_stats != false){
       c_run = ngsutils //without brackets

    } else if (params.with_fastqc != false && params.with_stats == false){
        c_run = fastqc
    } else {
        print ("Error: Please provide either --with_fastqc or --with_stats")
		System.exit(1)
    }
   prefetch | split | c_run 
}