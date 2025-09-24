nextflow.enable.dsl = 2 //which language do you want to use

process download_file {

//where should the output go?
publishDir "${projectDir}/output", mode: "copy", overwrite: true
//old version: publishDir "/home/lui/nextflow/output", mode: "copy", overwrite: true
//its also possible with "publishDir projectDir"


//what do you want to save permanently? hier brauche ich anführungszeichen für den path, sonst wird es falsch
output:
	path "batch1.fasta"

"""
	wget https://tinyurl.com/cqbatch1 -O batch1.fasta
"""
}

workflow {

	download_file()

}
