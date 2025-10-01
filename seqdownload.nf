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
4020500
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
    container "https://depot.galaxyproject.org/singularity/fastp%3A1.0.1--heae3180_0" //look for the corresponding container in galaxy
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

/*
fastp \
        --in1 ${inputvariable} \
        --out1 ${sample_id} \                               //generatats an output zip
        --cut_window_size ${params.cut_window_size} \
        --cut_mean_quality ${params.cut_mean_quality} \
        --length_required ${params.length_required} \
        --average_qual ${params.average_qual} \
        --html ${sample_id}_fastp.html \                    // generates html report
        --json ${sample_id}_fastp.json                         //generates json report


--in1 ${inputvariable}
Gibt die Eingabedatei mit den Sequenzreads an (FASTQ oder FASTQ.GZ).
Bei gepaarten Reads würdest du zusätzlich --in2 angeben.

--out1 ${sample_id}_trimmed.fastq.gz
Ausgabedatei für die gefilterten/qualitätsgeprüften Reads.
(Für Paired-End: auch --out2 notwendig.)


*/



workflow {
  
   accession = channel.fromPath(params.accession_number).splitText().map{it -> it.trim()}  
   //data_ch = prefetch(accession) | split

   prefetch_ch = prefetch(accession)
   split_ch = split(prefetch_ch) 

     if (params.with_fastp){
        fastp(split_ch)
    }
/*
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

  */


}


/*


dear chatgpt, please write a nextflow process that executes fastp on single end data. allow passing these params:

// --cut_window_size $params.something
// --cut_mean_quality $params.something
// --length_required $params.something
// --average_qual $params.something

ps: output html and json report


*/
