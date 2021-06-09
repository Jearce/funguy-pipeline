rule Rcorrector:
    input:
        r1 = "{species}/rnaseq/{species}_R1.fastq",
        r2 = "{species}/rnaseq/{species}_R2.fastq"
    output:
        "{species}/rnaseq/corrected/{species}_R1.cor.fq",
        "{species}/rnaseq/corrected/{species}_R2.cor.fq"
    threads: 20
    run:
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate rnaqc-env")
        shell("run_rcorrector.pl -1 {input.r1} -2 {input.r2} -od {wildcards.species}/rnaseq/corrected -t {threads}")
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda deactivate")

rule purge_unfixable:
    input:
        cor1 = "{species}/rnaseq/corrected/{species}_R1.cor.fq",
        cor2 = "{species}/rnaseq/corrected/{species}_R2.cor.fq"
    output:
        "{species}/rnaseq/corrected/unfixrm_{species}_R1.cor.fq",
        "{species}/rnaseq/corrected/unfixrm_{species}_R2.cor.fq"
    run:
        shell("python FilterUncorrectabledPEfastq.py -1 {input.cor1} -2 {input.cor2} -s {wildcards.species}")
        shell("mv unfixrm* {wildcards.species}/rnaseq/corrected")
        shell("mv rmunfix* {wildcards.species}/rnaseq/corrected")

rule trim_galore:
    input:
        purg1 = "{species}/rnaseq/corrected/unfixrm_{species}_R1.cor.fq",
        purg2 = "{species}/rnaseq/corrected/unfixrm_{species}_R2.cor.fq"
    params:
        qc_out="--outdir {species}/rnaseq/trimmed"
    output:
        "{species}/rnaseq/trimmed/{species}_val_1.fq",
        "{species}/rnaseq/trimmed/{species}_val_2.fq"
    run:
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate rnaqc-env")
        shell("trim_galore --paired --retain_unpaired --basename {wildcards.species} --fastqc_args {params.qc_out} "
        "--output_dir {species}/rnaseq/trimmed --length 36 -q 5 --cores 4 {input.purg1} {input.purg2}")
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda deactivate")

rule STAR:
    input:
        trim1 = "{species}/rnaseq/trimmed/{species}_val_1.fq",
        trim2 = "{species}/rnaseq/trimmed/{species}_val_2.fq"
    params:
        index_dir="{species}/rnaseq/mapped"
        genome_fasta="{species}/polished/best/polish.fasta"
        prefix="{species}/rnaseq/mapped/{species}."
    output:
        "{species}/rnaseq/mapped/{species}.Aligned.out.bam"
    threads: 24
    run:
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate star-env")
        shell("STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles "
        "{params.genome_fasta} --genomeSAindexNbases 12")
        shell("STAR --runThreadN {threads} --genomeDir {params.index_dir} --readFilesIn {input.trim1} {input.trim2} "
        "--outSAMtype BAM Unsorted SortedByCoordinate --outFileNamePrefix {params.prefix} --twopassMode Basic "
        "--outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 100000 --alignMatesGapMax 100000 "
        "--outFilterIntronMotifs RemoveNoncanonical --outReadsUnmapped Fastx")
