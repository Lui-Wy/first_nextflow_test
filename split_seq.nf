nextflow.enable.dsl = 2 
params.out = "${projectDir}/output"
params.temp = "${projectDir}/downloads" //generates a folder with downloads
params.prefix = "Seq_"

/*
for commenting out paragraphs
*/


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


process test_pass {
	//this means something like pass
	// I could also include: 
		// input: 
			//path input_variable 
		//output:
			//path "test_pass.txt"
	//that means that I could test something until this point in the working directory

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


process count_repeats {

	publishDir params.out, mode: "copy", overwrite: true

	input: 
		path input_variable 
	output:
		path "${input_variable.getSimpleName()}_count_repeats.txt" 
	"""
		grep -o "GCCGCG" ${input_variable}| wc -l > ${input_variable.getSimpleName()}_count_repeats.txt
	"""


//old: tail -n 1 ${input_variable} | wc -m > ${input_variable}_basecount.txt
}


process make_summary {
	publishDir params.out, mode: "copy", overwrite: true

	input: 
		path input_variable 
	output:
		path "summary.txt" 
	"""
		(for i in \$(ls ${input_variable} | sort); do echo \$i | cut -d "_" -f 1,2; cat \$i; done) > summary.txt
	"""
// ich muss das dollarzeichen im echo escapen! \$ , damit es im bash commando erhalten bleibt und nicht, dass hier groovy meint, das sei eine varibale
// alternative: for i in \$(ls ${input_variable}); do echo \$i; cat \$i; done > summary.txt
//sort \$i; (${input_variable} | sort)

/* 
(for i in \$(ls ${input_variable} | cut -d "_" -f 1,2 | sort); do echo \$i;  done) > summary.txt ----> works, but!
just works, if I dont need to use the cat function, because cat$i does not exist anymore after cut!
*/


}



workflow {

	
	download_output = download_file()
	count_seq (download_output)
	channel_split_seq = split_seq(download_output)

	channel_split_seq_flatten = channel_split_seq.flatten()

	count_bases(channel_split_seq_flatten)
	

	c_count_repeats = count_repeats(channel_split_seq_flatten)
	c_count_repeats.view()

	count_repeats_collect = c_count_repeats.collect()
	count_repeats_collect.view()

	make_summary(count_repeats_collect)
	

	// faster and more beautiful: 
	// download_file | split_seq | flatten | count_repeats | collect | make_summary


}
