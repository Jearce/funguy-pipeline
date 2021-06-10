#def get_cds_file_input(wildcards):
    #return config["cds"][wildcards.species]

rule ASGART:
    input:
        path = "{species}/polished/best/polish.fasta"
    output:
        "{species}_sd.json"
    shell:
        "asgart -k 101 --min-length 10000 --out {wildcards.species}_sd {input.path}"

rule ASGART_plot:
    input:
        path = "{species}_sd.json"
    output:
        "{species}/ASGART/{species}_sd.conf"
    run:
        shell("asgart-plot {input.path} --colorize by-position --out={wildcards.species}_sd circos")
        shell("mv *_sd* {wildcards.species}/ASGART")

rule EDTA:
    input:
        path = "{species}/polished/best/polish.fasta"                              #replace with purged final assembly
        #cds=get_cds_file_input
    output:
        "{species}/EDTA/polish.fasta.mod.EDTA.anno/polish.fasta.mod.EDTA.TEanno.out"
    params:
        cds = config["cds"]["{species}"]
    log:
        "{species}/EDTA/{species}_EDTA.log"
    threads: 20
    run:
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate edta && "
        "EDTA.pl --genome {input.path} --cds {params.cds} --sensitive 1 --anno 1 "
        "--evaluate 1 --species others -t {threads} 2> {log}")
        shell("mv polish.fasta* {wildcards.species}/EDTA")                       #modify with purged final assembly

rule Softmask: #Use *.mod.EDTA.anno/*.mod.EDTA.TEanno.out
    input:
        path = "{species}/polished/best/polish.fasta",                             #replace with purged final assembly
        lib = "{species}/EDTA/polish.fasta.mod.EDTA.anno/polish.fasta.mod.EDTA.TEanno.out"                      #modify with purged final assembly
    output:
        "{species}/polished/best/polish.fasta.masked"
    threads: 4
    run:
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate edta && "
        "make_masked.pl -genome {input.path} -minlen 80 -hardmask 0 -t {threads} -rmout {input.lib}")

rule BRAKER2: #Reference fungal protein sequences from OrthoDB
    input:
        path = "{species}/polished/best/polish.fasta.masked",
        proteins = "refseq/odb10_fungi_proteins.fasta",
        bam = "{species}/rnaseq/mapped/{species}.Aligned.out.bam"
    output:
        "{species}/BRAKER2/augustus.hints.aa"
    threads: 22
    run:
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate braker && "
        "braker.pl --species={wildcards.species} --genome={input.path} --prot_seq={input.proteins} --softmasking "
        "--etpmode --fungus --cores {threads} --workingdir={wildcards.species}/BRAKER2")
