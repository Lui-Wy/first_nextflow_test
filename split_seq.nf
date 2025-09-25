nextflow.enable.dsl = 2 
params.out = "${projectDir}/output"
params.temp = "${projectDir}/downloads" //generates a folder with downloads
params.prefix = "Seq_"

process download_file {
	publishDir params.out, mode: "copy", overwrite: true
	storeDir params.temp //looks for the name wich I give under wget -O in the download folder, if its there, download is skipped
	output:
		path "batch1.fasta" 

	"""
		wget https://tinyurl.com/cqbatch1 -O batch1.fasta
	"""
}


process count_seq{
	publishDir params.out, mode: "copy", overwrite: true

	input: 
		path input_variable 
	output:
		path "count_seq.txt"
	"""
		grep ">" ${input_variable} | wc -l > count_seq.txt
	"""
		
}


process test {
	"""
	
	"""
		
}


process split_seq {
	publishDir params.out, mode: "copy", overwrite: true

	input: 
		path input_variable 
	output:
		path "${params.prefix}*.fasta"  //has to be the name, which comes from the shell (output of(split -d -l 2 --additional-suffix .fasta ${output_variable} Seq_ ))

	"""
		split -d -l 2 --additional-suffix .fasta ${input_variable} ${params.prefix}
	"""
	// split generates files, so I dont need the addition liek > "*.txt"
}


process count_bases {

	publishDir params.out, mode: "copy", overwrite: true

	input: 
		path input_variable 
	output:
		path "${input_variable.getSimpleName()}_basecount.txt" //getSimpleName is a method for paths(file)
	"""
		tail -n 1 ${input_variable} | wc -m > ${input_variable.getSimpleName()}_basecount.txt 
	"""


//old: tail -n 1 ${input_variable} | wc -m > ${input_variable}_basecount.txt
}





workflow {

	download_output = download_file()
	count_seq (download_output)
	channel_split_seq = split_seq(download_output)
	channel_split_seq_flatten = channel_split_seq.flatten()
	count_bases(channel_split_seq_flatten)


}
