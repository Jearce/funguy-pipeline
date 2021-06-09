rule ASGART:
    input:
        path="{species}/polished/best/polish.fasta"
    output:
        "{species}_sd.json"
    shell:
        "asgart -k 101 --min-length 10000 --out {wildcards.species}_sd {input.path}"

rule ASGART_plot:
    input:
        path="{species}_sd.json"
    output:
        "{species}/ASGART/{species}_sd.conf"
    run:
        shell("asgart-plot {input.path} --colorize by-position --out={wildcards.species}_sd circos")
        shell("mv *_sd* {wildcards.species}/ASGART")

rule EDTA:
    input:
        path="{species}/polished/best/polish.fasta",                             #replace with purged final assembly
        #cds=get_cds_file_input
    output:
        "{species}/EDTA/polish.fasta.mod.EDTA.TElib.fa"
    params:
        cds=config["cds"]["{species}"]
    log:
        "{species}/{species}_EDTA.log"
    threads: 20
    run:
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate EDTA-env")
        shell("EDTA.pl --genome {input.path} --cds {params.cds} --sensitive 1 --anno 1 "
        "--evaluate 1 --species others -t {threads} 2> {log}")
        shell("mv polish.fasta* {wildcards.species}/EDTA")                       #modify with purged final assembly
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda deactivate")

rule RepeatMasker: #Use repeat library constructed from EDTA
    input:
        path="{species}/polished/best/polish.fasta",                             #replace with purged final assembly
        lib="{species}/EDTA/polish.fasta.mod.EDTA.TElib.fa"                      #modify with purged final assembly
    output:
        "{species}/softmasked/polish.fasta.masked"
    threads: 4
    run:
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate EDTA-env")
        shell("RepeatMasker -s -xsmall -engine ncbi -lib {input.lib} "
        "-pa {threads} -dir {wildcards.species}/softmasked {input.path}")
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda deactivate")

rule BRAKER2: #Reference fungal protein sequences from OrthoDB
    input:
        path="{species}/softmasked/polish.fasta.masked",
        proteins="refseq/odb10_fungi_proteins.fasta"
    output:
        "{species}/BRAKER2/braker.gtf"
    threads: 22
    run:
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda activate braker-env")
        shell("braker.pl --species={wildcards.species} --genome={input.path} "
        "--prot_seq={input.proteins} --softmasking --epmode --fungus "
        "--cores {threads} --workingdir={wildcards.species}/BRAKER2")
        shell("set +eu && PS1=dummy && . $(conda info --base)/etc/profile.d/conda.sh && conda deactivate")
