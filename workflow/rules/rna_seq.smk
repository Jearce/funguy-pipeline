rule Rcorrector:
    input:
        r1 = "{species}/rnaseq/{species}_R1.fastq",
        r2 = "{species}/rnaseq/{species}_R2.fastq"
    output:
        temp("{species}/rnaseq/corrected/{species}_R1.cor.fq"),
        temp("{species}/rnaseq/corrected/{species}_R2.cor.fq")
    threads: 20
    conda:
        "envs/rcorrector.yaml"
    shell:
        "run_rcorrector.pl -1 {input.r1} -2 {input.r2} -od {wildcards.species}/rnaseq/corrected -t {threads}"

rule purge_unfixable:
    input:
        cor1 = "{species}/rnaseq/corrected/{species}_R1.cor.fq",
        cor2 = "{species}/rnaseq/corrected/{species}_R2.cor.fq"
    output:
        temp("{species}/rnaseq/corrected/unfixrm_{species}_R1.cor.fq"),
        temp("{species}/rnaseq/corrected/unfixrm_{species}_R2.cor.fq")
    run:
        shell("python FilterUncorrectabledPEfastq.py -1 {input.cor1} -2 {input.cor2} -s {wildcards.species}")
        shell("mv unfixrm* {wildcards.species}/rnaseq/corrected")
        shell("mv rmunfix* {wildcards.species}/rnaseq/corrected")

rule trim_galore:
    input:
        purg1 = "{species}/rnaseq/corrected/unfixrm_{species}_R1.cor.fq",
        purg2 = "{species}/rnaseq/corrected/unfixrm_{species}_R2.cor.fq"
    params:
        qc_out = "--outdir {species}/rnaseq/trimmed"
    output:
        "{species}/rnaseq/trimmed/{species}_val_1.fq",
        "{species}/rnaseq/trimmed/{species}_val_2.fq"
    conda:
        "envs/trim_galore.yaml"
    shell:
        "trim_galore --paired --retain_unpaired --basename {wildcards.species} --output_dir {wildcards.species}/rnaseq/trimmed "
        "--stringency 3 --length 36 -q 20 --cores 4 --fastqc_args {params.qc_out} {input.purg1} {input.purg2}"

rule STAR:
    input:
        trim1 = "{species}/rnaseq/trimmed/{species}_val_1.fq",
        trim2 = "{species}/rnaseq/trimmed/{species}_val_2.fq",
        genome_fasta = "{species}/polished/best/polish.fasta"
    params:
        index_dir = "{species}/rnaseq/mapped",
        prefix = "{species}/rnaseq/mapped/{species}."
    output:
        "{species}/rnaseq/mapped/{species}.Aligned.out.bam"
    conda:
        "envs/star.yaml"
    threads: 24
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles "
        "{input.genome_fasta} --genomeSAindexNbases 12 && "
        "STAR --runThreadN {threads} --genomeDir {params.index_dir} --readFilesIn {input.trim1} {input.trim2} "
        "--outSAMtype BAM Unsorted SortedByCoordinate --outFileNamePrefix {params.prefix} --twopassMode Basic "
        "--outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 10000 --alignMatesGapMax 1000000 "
        "--outFilterIntronMotifs RemoveNoncanonical --outReadsUnmapped Fastx"
