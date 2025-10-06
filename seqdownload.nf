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
    storeDir params.temp

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
    
    publishDir params.out, mode: "copy", overwrite: true
    input:
        path inputvariable
    output:
        path "${sample_id}_trimmed.fastq", emit: trimmed // - multiqc just needs the json file, so I give here aliases
        path "${sample_id}_fastp.html", emit: html 
        path "${sample_id}_fastp.json", emit: json //alternative: I could use aliases with ", emit: json"
    script:
        sample_id = inputvariable.getSimpleName()

        """
        fastp \
        --in1 ${inputvariable} \
        --out1 ${sample_id}_trimmed.fastq  \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json 
        """
}



process multi_qc {
 container "https://depot.galaxyproject.org/singularity/multiqc%3A1.9--pyh9f0ad1d_0"
   
    publishDir params.out, mode: "copy", overwrite: true
    input:
        path inputvariable
    output:
        path "multiqc_*" // not "${sample_id}.html" --> has another name, it always puts out
   
        """
        multiqc .
        """


    }


// von Max 

process spades {
    publishDir params.out, mode: 'copy', overwrite: true
    container "https://depot.galaxyproject.org/singularity/spades%3A4.2.0--h8d6e82b_2"
    input:
        path fastqFiles
    output:
        path "${fastqFiles.getSimpleName()}"
     script:
        """
            spades.py -s ${fastqFiles} -o ${fastqFiles.getSimpleName()}
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
        fastp_ch = fastp(split_ch)
    } else {
        print ("Attention: If you want to do a fastp- analysis, please provide --with_fastp")
        //System.exit(1)                                                            //bedeutet auch, sobald ich --fastp nicht setze, bricht der run ab! (auch kein fastqc/stats)
    }


    multi_qc(fastp_ch.json.collect()) //needs to collect all three files

    spades (fastp_ch.trimmed)



}

