import gzip
from Bio import SeqIO

configfile: "../config/assembly/config.yaml"

def get_genome_size(wildcards):
	return config[wildcards.species]["genomeSize"]
	
def fastq_count(seqfile):
	allbp=0
	with gzip.open(seqfile, "rt") as handle:
		for line in SeqIO.parse(handle, 'fastq'):
			allbp+=len(line.seq)
	return allbp

rule canu_correction: #produce the corrected reads in {species}/corrected_reads dir then mv it to {species} dir and rename to {output}
	input:
		path="{species}/nanopore/{species}.fastq.gz"
	output:
		"{species}/corrected_{species}_nano.fasta.gz"
	params:
		gsize=get_genome_size
	threads: 30
	run:
		if (fastq_count(input.path)/(float(params.gsize)*(10**6)) <= 20):
			print('LOW COVERAGE AT CORRECTION STEP. USING ALL READS IN NANOPORE INPUTS')
			shell("cat {input.path} > {output}")
		else:
			shell("canu -correct -p canu genomeSize={params.gsize}m stopOnLowCoverage=1 "
			"-d {wildcards.species}/corrected_reads -nanopore {input.path} -maxthreads={threads} -maxmemory=60 -useGrid=0")
			shell("mv {wildcards.species}/corrected_reads/canu.correctedReads.fasta.gz {output}")

rule canu_trim: #produce the trimmed reads in {species}/trimmed_reads dir then mv it to {species} dir and rename to {output}
	input:
		path="{species}/corrected_{species}_nano.fasta.gz",
		raw="{species}/nanopore/{species}.fastq.gz"
	output:
		"{species}/trimmed_corr_{species}_nano.fasta.gz"
	params:
		gsize=get_genome_size
	threads: 30
	run:
		if (fastq_count(input.raw)/(float(params.gsize)*(10**6)) <= 20):
			print('LOW COVERAGE AT TRIMMING STEP. USING ALL READS IN NANOPORE INPUTS')
			shell("cat {input.path} > {output}")
		else:
			shell("canu -trim -p canu genomeSize={params.gsize}m stopOnLowCoverage=1 "
			"-corrected -d {wildcards.species}/trimmed_reads -nanopore {input.path} -maxthreads={threads} -maxmemory=60 -useGrid=0")
			shell("mv {wildcards.species}/trimmed_reads/canu.trimmedReads.fasta.gz {output}")

rule canu_assemble:	#requires trimmed and corrected reads from canu. Output to {species}/canu_out.
					#Assembled genome is moved to {species}/drafts folder
	input:
		path="{species}/trimmed_corr_{species}_nano.fasta.gz",
		raw="{species}/nanopore/{species}.fastq.gz"
	output:
		"{species}/drafts/canu_{species}_contigs.fasta"
	params:
		gsize=get_genome_size
	threads: 20
	run:
		shell("if [ -d {wildcards.species}/drafts ]; then echo drafts folder already exists; else mkdir {wildcards.species}/drafts; fi")
		if (fastq_count(input.raw)/(float(params.gsize)*(10**6)) <= 20):
			shell("canu -p {wildcards.species} genomeSize={params.gsize}m corMinCoverage=0 correctedErrorRate=0.25 -maxthreads={threads} -maxmemory=60 "
			"corMhapSensitivity=high minReadLength=600 -d {wildcards.species}/canu_out -nanopore {input.path} -useGrid=0 "
			"&& cp {wildcards.species}/canu_out/{wildcards.species}.contigs.fasta {wildcards.species}/drafts/")
		else:
			shell("canu -p {wildcards.species} genomeSize={params.gsize}m stopOnLowCoverage=1 -maxthreads={threads} -maxmemory=60 "
			"-trimmed -corrected -d {wildcards.species}/canu_out -nanopore {input.path} -useGrid=0 "
			"&& cp {wildcards.species}/canu_out/{wildcards.species}.contigs.fasta {wildcards.species}/drafts/")
		shell("mv {wildcards.species}/drafts/{wildcards.species}.contigs.fasta {output}")
		
rule flye_assemble: #requires trimmed and corrected reads from canu. Output to flye_out. 
					#Assembled genome is moved to {species}/drafts folder
	input:
		path="{species}/trimmed_corr_{species}_nano.fasta.gz",
		raw="{species}/nanopore/{species}.fastq.gz"
	output:
		"{species}/drafts/flye_{species}_assembly.fasta"
	threads: 15
	params:
		gsize=get_genome_size
	run:
		shell("if [ -d {wildcards.species}/drafts ]; then echo drafts folder already exists; else mkdir {wildcards.species}/drafts; fi")
		if (fastq_count(input.raw)/(float(params.gsize)*(10**6)) <= 20):
			print('LOW COVERAGE FLYE INPUT. USING RAW DATA')
			shell("flye -g {params.gsize}m -t {threads} "
			"-o {wildcards.species}/flye_out --nano-raw {input.raw}"
			"&& cp {wildcards.species}/flye_out/assembly.fasta {wildcards.species}/")
		else:
			shell("flye -g {params.gsize}m -t {threads} "
			"-o {wildcards.species}/flye_out --nano-corr {input.path}"
			"&& cp {wildcards.species}/flye_out/assembly.fasta {wildcards.species}/")
		shell("mv {wildcards.species}/assembly.fasta {output}")
		
rule zipping_short_1:
	input:
		short1="{species}/illumina/{species}_R1.fastq"
	output:
		out1="{species}/illumina/{species}_R1.fastq.gz","
	run:
		print('Zipping first short read')
		shell("if [ ! -f {input.short1}.gz ]; then gzip -c {input.short1} > {output.out1}; fi")
		
rule zipping_short_2:
	input:
		short2="{species}/illumina/{species}_R2.fastq"
	output:
		out2="{species}/illumina/{species}_R2.fastq.gz","
	run:
		print('Zipping first short read')
		shell("if [ ! -f {input.short2}.gz ]; then gzip -c {input.short2} > {output.out2}; fi")
		
rule zipping_long:
	input:
		long="{species}/nanopore/{species}.fastq"
	output:
		out3="{species}/nanopore/{species}.fastq.gz","
	run:
		print('Zipping long read')
		shell("if [ ! -f {input.long}.gz ]; then gzip -c {input.long} > {output.out3}; fi")

rule wengan_assemble: 	#produce all files in cwd. After assembly, move everything to wengan_out
						#assembled genome is moved to {species} folder
						#input MUST be .gz Other file types may result in errors
	input:
		nano="{species}/nanopore/{species}.fastq.gz",
		short1="{species}/illumina/{species}_R1.fastq.gz",
		short2="{species}/illumina/{species}_R2.fastq.gz"
	output:
		"{species}/drafts/wengan_{species}_assembly.fasta"
	threads: 5
	params:
		gsize=get_genome_size
	shadow: "full"
	run:
		shell("if [ -d {wildcards.species}/drafts ]; then echo drafts folder already exists; else mkdir {wildcards.species}/drafts; fi")
		shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate wengan-runtime-env "
		"&& perl $WG -x ontraw -a M -s {input.short1},{input.short2} "
		"-l {input.nano} -p wengan_{wildcards.species} -t {threads} -g {params.gsize} "
		"&& mkdir {wildcards.species}/wengan_out/")
		shell("mv wengan_{wildcards.species}* {wildcards.species}/wengan_out/")
		shell("mv {wildcards.species}/wengan_out/wengan_{wildcards.species}.SPolished.asm.wengan.fasta {output}")
