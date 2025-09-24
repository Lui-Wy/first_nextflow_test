nextflow.enable.dsl = 2 //which language do you want to use

process download_file {
	//where should the output go?
	publishDir "${projectDir}/output", mode: "copy", overwrite: true
	//old version: publishDir "/home/lui/nextflow/output", mode: "copy", overwrite: true
	//its also possible with "publishDir projectDir"


	//what do you want to save permanently? hier brauche ich anführungszeichen für den path, sonst wird es falsch
	output:
		path "batch1.fasta" //path means here something like file

	"""
		wget https://tinyurl.com/cqbatch1 -O batch1.fasta
	"""
}

process count_seq {
	publishDir "${projectDir}/output", mode: "copy", overwrite: true

	input: 
		path output_variable // define that the output of 1. process (download) is now input for this process (in combination with the workflow download_output = download_file()
	                         //count_seq(download_output))
	output:
		path "number_of_seq.txt"

	"""
		grep ">" ${output_variable} | wc -l > number_of_seq.txt
	"""
}


workflow {

	download_output = download_file()//I need to put the output into a new variable, because I can run a process just once
	//if I need the output from download once again in a second process (not yet here) (same for the use in a pipe)
	count_seq(download_output)

}
