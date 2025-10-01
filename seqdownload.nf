nextflow.enable.dsl = 2 

params.temp = "${projectDir}/downloads" 
params.out = "${projectDir}/output"
params.with_fastqc = false 
params.with_stats = false

params.accession_number = "${projectDir}/accessions.txt" // entspricht: "/path/to/accessions.txt"



process prefetch {
    //storeDir params.temp

    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"

    input:
      val accession
    output:
        path "${accession}/${accession}.sra" 
        // path "${accession_number}/${accession_number}.sra" would mean:
        // path/to/accessions.txt//path/to/accessions.txt.sra" /
      

    """
    prefetch ${accession}
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

}


process ngsutils {
    container "https://depot.galaxyproject.org/singularity/ngsutils%3A0.5.9--py27heb79e2c_4"
    //storeDir params.temp 
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
        path "${inputvariable.getSimpleName()}_fastqc.zip"
        path "${inputvariable.getSimpleName()}_fastqc.html"

    """
    fastqc -o . -f fastq "${inputvariable}" 
    """
}


workflow {
  
   accession = channel.fromPath(params.accession_number).splitText().map{it -> it.trim()}
   prefetch_ch = prefetch(accession)
   split_ch = split(prefetch_ch) 

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


/*

SRR1777174
SRR12718173
SRR12426925

accession = channel.fromPath(params.accessions).splitText().map{it -> it.trim()}


dear chatgpt, please write a nextflow process that executes fastp on single end data. allow passing these params:

// --cut_window_size $params.something
// --cut_mean_quality $params.something
// --length_required $params.something
// --average_qual $params.something

ps: output html and json report


*/
