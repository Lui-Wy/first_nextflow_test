nextflow.enable.dsl = 2 

params.temp = "${projectDir}/downloads" 
params.out = "${projectDir}/output"

//params for conditionals at the end:
params.with_fastqc = false 
params.with_stats = false

//params for multiple input accession numbers:
params.accession_number = "${projectDir}/accessions.txt" // entspricht: "/path/to/accessions.txt"

//params for fastp-trimming:
params.cut_window_size = 4 // set default values from the man fastp page 
params.cut_mean_quality = 20
params.length_required = 50
params.average_qual = 0

//making fastq optional
params.with_fastp = false


process prefetch {
    //storeDir params.temp

    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"

    input:
      val accession
    output:
        path "${accession}/${accession}.sra" 
            
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



process fastp {
    container "https://depot.galaxyproject.org/singularity/fastp%3A1.0.1--heae3180_0" 
    //storeDir params.temp
    publishDir params.out, mode: "copy", overwrite: true
    input:
        path inputvariable
    output:
        path "${sample_id}_trimmed.fastq.gz"
        path "${sample_id}_fastp.html" 
        path "${sample_id}_fastp.json"
    script:
        sample_id = inputvariable.getSimpleName()

        """
        fastp \
        --in1 ${inputvariable} \
        --out1 ${sample_id}_trimmed.fastq.gz  \
        --cut_window_size ${params.cut_window_size} \
        --cut_mean_quality ${params.cut_mean_quality} \
        --length_required ${params.length_required} \
        --average_qual ${params.average_qual} \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json 
        """
}




workflow {
  
   accession = Channel.fromPath(params.accession_number).splitText().map{it -> it.trim()}  
   //data_ch = prefetch(accession) | split

   prefetch_ch = prefetch(accession)
   split_ch = split(prefetch_ch) 

    if (params.with_stats){
       ngsutils(split_ch)
    } 
    if (params.with_fastqc){
        fastqc(split_ch)
    } 
    if (params.with_fastqc == false && params.with_stats == false) {
        print ("Attention for fastqc or stats, please provide either --with_fastqc or --with_stats")
		//System.exit(1)                                                              //wenn ich das hier schreibe, wird der folgende if block gar nicht ausgeführt!
    }

    if (params.with_fastp){
        fastp(split_ch)
    } else {
        print ("Attention: If you want to do a fastp- analysis, please provide either --with_fastp")
        //System.exit(1)                                                            //bedeutet auch, sobald ich --fastp nicht setze, bricht der run ab! (auch kein fastqc/stats)
    }


    /*
    if (params.with_stats){
       ngsutils (split_ch)
    } else if (params.with_fastqc){
        fastqc (split_ch)
    } else if (params.with_fastqc && params.with_stats){ //wenn ich das so schreibe, wird diese bedinnung nie erfüllt, da ich das ja schon so im ersten if (params.with_stats) habe
        ngsutils (split_ch)
        fastqc (split_ch)
    } else {
        print ("Error: Please provide either --with_fastqc or --with_stats")
		System.exit(1)
    }

  */


}

