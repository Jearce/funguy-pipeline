configfile: "../../config/assembly/config.yaml"

def get_genome_size(wildcards):
	return config[wildcards.species]["genomeSize"]

rule canu_correction: #produce the corrected reads in {species}/corrected_reads dir then mv it to {species} dir and rename to {output}
	input:
		path="{species}/nanopore/{species}.fastq"
	output:
		"{species}/corrected_{species}_nano.fasta.gz"
	params:
		gsize=get_genome_size
	run:
		shell("canu -correct -p canu genomeSize={params.gsize}m "
		"-d {wildcards.species}/corrected_reads -nanopore {input.path} "
		"&& mv {wildcards.species}/corrected_reads/canu.correctedReads.fasta.gz {output}")

rule canu_trim: #produce the trimmed reads in {species}/trimmed_reads dir then mv it to {species} dir and rename to {output}
	input:
		path="{species}/corrected_{species}_nano.fasta.gz"
	output:
		"{species}/trimmed_corr_{species}_nano.fasta.gz"
	params:
		gsize=get_genome_size
	run:
		shell("canu -trim -p canu genomeSize={params.gsize}m "
		"-corrected -d {wildcards.species}/trimmed_reads -nanopore {input.path} "
		"&& mv {wildcards.species}/trimmed_reads/canu.trimmedReads.fasta.gz {output}")

rule canu_assemble:	#requires trimmed and corrected reads from canu. Output to {species}/canu_out.
					#Assembled genome is moved to {species}/drafts folder
	input:
		path="{species}/trimmed_corr_{species}_nano.fasta.gz"
	output:
		"{species}/drafts/canu_{species}_contigs.fastq"
	params:
		gsize=get_genome_size
	threads: 15
	run:
		shell("if [ -d {wildcards.species}/drafts ]; then echo drafts folder already exists; else mkdir {wildcards.species}/drafts; fi")
		shell("canu -p {wildcards.species} genomeSize={params.gsize}m -maxthreads={threads} "
		"-trimmed -corrected -d {wildcards.species}/canu_out -nanopore {input.path}"
		"&& cp {wildcards.species}/canu_out/{wildcards.species}.contigs.fasta {wildcards.species}/drafts/")
		shell("mv {wildcards.species}/drafts/{wildcards.species}.contigs.fasta {output}")
		
rule flye_assemble: #requires trimmed and corrected reads from canu. Output to flye_out. 
					#Assembled genome is moved to {species}/drafts folder
	input:
		path="{species}/trimmed_corr_{species}_nano.fasta.gz"
	output:
		"{species}/drafts/flye_{species}_assembly.fasta"
	threads: 10
	params:
		gsize=get_genome_size
	run:
		shell("if [ -d {wildcards.species}/drafts ]; then echo drafts folder already exists; else mkdir {wildcards.species}/drafts; fi")
		shell("flye -g {params.gsize}m -t {threads} "
		"-o {wildcards.species}/flye_out --nano-corr {input.path}"
		"&& cp {wildcards.species}/flye_out/assembly.fasta {wildcards.species}/")
		shell("mv {wildcards.species}/assembly.fasta {output}")
		
rule zipping_file:
	input:
		short1="{species}/illumina/{species}_R1.fastq",
		short2="{species}/illumina/{species}_R2.fastq"
	output:
		out1="{species}/illumina/{species}_R1.fastq.gz",
		out2="{species}/illumina/{species}_R2.fastq.gz"
	run:
		shell("gzip {input.short1} && gzip {input.short2}")

rule wengan_assemble: 	#produce all files in cwd. After assembly, move everything to wengan_out
						#assembled genome is moved to {species} folder
						#input MUST be .gz Other file types may result in errors
						#Seem to require high coverage ~30X minimum
	input:
		nano="{species}/corrected_{species}_nano.fasta.gz",
		short1="{species}/illumina/{species}_R1.fastq.gz",
		short2="{species}/illumina/{species}_R2.fastq.gz"
	output:
		"{species}/drafts/wengan_{species}_assembly.fasta"
	threads: 5
	params:
		gsize=get_genome_size
	run:
		shell("if [ -d {wildcards.species}/drafts ]; then echo drafts folder already exists; else mkdir {wildcards.species}/drafts; fi")
		shell("wengan.pl -x ontraw -a M -s {input.short1},{input.short2} "
		"-l {input.nano} -p wengan_{wildcards.species} -t {threads} -g {gsize} "
		"&& mkdir {wildcards.species}/wengan_out/")
		shell("mv wengan_{wildcards.species}* {wildcards.species}/wengan_out/")
		shell("mv {wildcards.species}/wengan_out/wengan_{wildcards.species}.SPolished.asm.wengan.fasta {output}")
