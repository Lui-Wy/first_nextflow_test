nextflow.enable.dsl = 2 

params.temp = "${projectDir}/downloads" 
params.out = "${projectDir}/output"



process prefetch {
    storeDir params.temp
    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"

    output:
        path "SRR1777174" //alternative just the file in the folder so: "SRR1777174/SRR1777174.sra"

"""
prefetch SRR1777174
"""
}

process split {
    storeDir params.temp
    publishDir params.out, mode: "copy", overwrite: true

    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"
    input:
        path inputvariable
    output:
        path "*.fastq"

"""
fastq-dump --split-3 "${inputvariable}" 
"""
//split-3 does already generate an output - so I don´t need a redirection like > "test.txt"

}


process ngsutils {
    container "https://depot.galaxyproject.org/singularity/ngsutils%3A0.5.9--py27heb79e2c_4"
    storeDir params.temp
    publishDir params.out, mode: "copy", overwrite: true
   input:
        path inputvariable
    output:
        path "stats.txt"

"""
fastqutils stats "${inputvariable}" > "stats.txt"
"""

}

process fastqc {
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    storeDir params.temp
    publishDir params.out, mode: "copy", overwrite: true
   input:
        path inputvariable
    output:
        path "*"

"""
fastqc -o . -f fastq SRR1777174.fastq 
"""

}
//    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] [-c contaminant file] seqfile1 .. seqfileN



workflow {
    c_prefetch = prefetch()
    c_split = split(c_prefetch)
    ngsutils(c_split)

}