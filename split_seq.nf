nextflow.enable.dsl = 2 
params.out = "${projectDir}/output"
params.temp = "${projectDir}/downloads" //generates a folder with downloads

process download_file {
	publishDir params.out, mode: "copy", overwrite: true
	storeDir params.temp //looks for the name wich I give under wget -O in the download folder, if its there, download is skipped
	output:
		path "batch1.fasta" 

	"""
		wget https://tinyurl.com/cqbatch1 -O batch1.fasta
	"""
}

process split_seq {
	publishDir params.out, mode: "copy", overwrite: true

	input: 
		path output_variable 
	output:
		path "Seq_*.fasta"  //has to be the name, which comes from the shell (output of(split -d -l 2 --additional-suffix .fasta ${output_variable} Seq_ ))

	"""
		split -d -l 2 --additional-suffix .fasta ${output_variable} Seq_  
	"""
	//brauch hier nicht > "*.txt"
}


workflow {

	download_output = download_file()
	split_seq(download_output)

}
