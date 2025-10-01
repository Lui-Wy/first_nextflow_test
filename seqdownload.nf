nextflow.enable.dsl = 2 

params.temp = "${projectDir}/downloads" 
params.out = "${projectDir}/output"
params.with_fastqc = false 
params.with_stats = false

params.accession_number = "SRR1777174"


process prefetch {
    //storeDir params.temp

    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"

    //input:
      //  path 
    output:
        //path "${params.accession_number}/${params.accession_number}.sra" //does not work - take another look here - maybe now
        path "${params.accession_number}" 

    """
    prefetch ${params.accession_number}
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
        //path "${inputvariable.getSimpleName()}*" // I generate two output files, so I need the "*" here, but also puts out the input

        path "${inputvariable.getSimpleName()}_fastqc.zip"
        path "${inputvariable.getSimpleName()}_fastqc.html"

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

    /*
    if (params.with_fastqc == false && params.with_stats != false){
       c_run = ngsutils //without brackets //run_ch

    } else if (params.with_fastqc != false && params.with_stats == false){
        c_run = fastqc
   
    } else {
        print ("Error: Please provide either --with_fastqc or --with_stats")
		System.exit(1)
    }
   */

    //prefetch  | split | c_run



    // this part here is because the pipe cant deal with, if both params --with_fastqc and --with_stats are given
   prefetch_ch = prefetch()
   split_ch = split(prefetch_ch) 

   /* old version, but two long
    if (params.with_fastqc == false && params.with_stats != false){
       ngsutils (split_ch)

    } else if (params.with_fastqc != false && params.with_stats == false){
        fastqc (split_ch)

    } else if (params.with_fastqc != false && params.with_stats != false){
        ngsutils (split_ch)
        fastqc (split_ch)
    } else {
        print ("Error: Please provide either --with_fastqc or --with_stats")
		System.exit(1)
    }
   */


    //still not quit right?
    if (params.with_stats){
       ngsutils (split_ch)
    } else if (params.with_fastqc){
        fastqc (split_ch)
    } else if (params.with_fastqc && params.with_stats){
        ngsutils (split_ch)
        fastqc (split_ch)
    } else {
        print ("Error: Please provide either --with_fastqc or --with_stats")
		System.exit(1)
    }



}