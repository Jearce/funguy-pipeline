##Genome size is stored in est_genomesize.txt as a number divided by 1million.
##E.g 5m bp is 5, 1g bp is 1000

#TODO: rename the prot_input for eggnog | Check if eggnog work from start to finish by itself
#TODO: Edit the threads CPU for rules

##Input format for assembly: nanopore_{species}.fasta	trimmed_{species}_illumina1.fastq.gz	trimmed_{species}_illumina2.fastq.gz

def read_file(species_id,filename):
	#with open(species_id+'/'+filename,"r") as file:
	#	gs=file.readline().strip()
	return open(species_id+'/'+filename,"r").readline().strip()

rule canu_correction: #produce the corrected reads in {species}/corrected_reads dir then mv it to {species} dir and rename to {output}
	input:
		path="{species}/nanopore_{species}.fasta" ##Need to change to fastq if input is different
	output:
		"{species}/corrected_{species}_nano.fasta.gz"
#	params:
#		gsize=read_genome_size
	run:
		gsize=read_file(wildcards.species,'est_genomesize.txt')
		shell("canu -correct -p canu genomeSize={gsize}m "
		"-d {wildcards.species}/corrected_reads -nanopore {input.path}"
		"&& mv {wildcards.species}/corrected_reads/canu.correctedReads.fasta.gz {output}")
		#shell("touch {output}")
		
rule canu_trim: #produce the trimmed reads in {species}/trimmed_reads dir then mv it to {species} dir and rename to {output}
	input:
		path="{species}/corrected_{species}_nano.fasta.gz"
	output:
		"{species}/trimmed_corr_{species}_nano.fasta.gz"
	run:
		gsize=read_file(wildcards.species,'est_genomesize.txt')
		shell("canu -trim -p canu genomeSize={gsize}m "
		"-corrected -d {wildcards.species}/trimmed_reads -nanopore {input.path}"
		"&& mv {wildcards.species}/trimmed_reads/canu.trimmedReads.fasta.gz {output}")
		#shell("touch {output}")
		
rule canu_assemble:	#requires trimmed and corrected reads from canu. Output to flye_out.
					#Assembled genome is moved to {species} folder
	input:
		path="{species}/trimmed_corr_{species}_nano.fasta.gz"
	output:
		"{species}/drafts/canu_{species}_contigs.fasta"
	run:
		shell("if [ -d {wildcards.species}/drafts ]; then echo drafts folder already exists; else mkdir {wildcards.species}/drafts; fi")
		gsize=read_file(wildcards.species,'est_genomesize.txt')
		shell("canu -p {wildcards.species} genomeSize={gsize}m "
		"-trimmed -corrected -d {wildcards.species}/canu_out -nanopore {input.path}"
		"&& cp {wildcards.species}/canu_out/{wildcards.species}.contigs.fasta {wildcards.species}/drafts/")
		shell("mv {wildcards.species}/drafts/{wildcards.species}.contigs.fasta {output}")
		
rule flye_assemble: #requires trimmed and corrected reads from canu. Output to flye_out. 
					#Assembled genome is moved to {species} folder
	input:
		path="{species}/trimmed_corr_{species}_nano.fasta.gz"
	output:
		"{species}/drafts/flye_{species}_assembly.fasta"
	threads: 4
	run:
		shell("if [ -d {wildcards.species}/drafts ]; then echo drafts folder already exists; else mkdir {wildcards.species}/drafts; fi")
		gsize=read_file(wildcards.species,'est_genomesize.txt')
		shell("flye -g {gsize}m -t {threads} "
		"-o {wildcards.species}/flye_out --nano-corr {input.path}"
		"&& cp {wildcards.species}/flye_out/assembly.fasta {wildcards.species}/")
		shell("mv {wildcards.species}/assembly.fasta {output}")
		#shell("touch {output}")
		
rule wengan_assemble: 	#produce all files in cwd. After assembly, move everything to wengan_out
						#assembled genome is moved to {species} folder
						#input MUST be .gz Other file types may result in errors
						#Seem to require high coverage ~30X minimum
	input:
		nano="{species}/corrected_{species}_nano.fasta.gz",
		short1="{species}/trimmed_{species}_illumina1.fastq.gz",
		short2="{species}/trimmed_{species}_illumina2.fastq.gz"
	output:
		"{species}/drafts/wengan_{species}_assembly.fasta"
	threads: 4
	run:
		shell("if [ -d {wildcards.species}/drafts ]; then echo drafts folder already exists; else mkdir {wildcards.species}/drafts; fi")
		gsize=read_file(wildcards.species,'est_genomesize.txt')
		shell("perl $WG -x ontraw -a M -s {input.short1},{input.short2} "
		"-l {input.nano} -p wengan_{wildcards.species} -t {threads} -g {gsize} "
		"&& mkdir {wildcards.species}/wengan_out/")
		shell("mv wengan_{wildcards.species}* {wildcards.species}/wengan_out/")
		shell("mv {wildcards.species}/wengan_out/wengan_{wildcards.species}.SPolished.asm.wengan.fasta {output}")
		
rule test_all:
	input:
		"{species}/wengan_{species}_assembly.fasta",
		"{species}/flye_{species}_assembly.fasta",
		"{species}/canu_{species}_contigs.fasta"
	output:
		"{species}/complete.txt"
	shell:
		"cat {input} > {output}"
