import subprocess

from pathlib import Path

configfile: "../config/assembly/config.yaml"

def get_lineage(wildcards):
  return config[wildcards.species]["lineage"]

rule reciprocal_polishing:
  input:
    r1 = "{species}/illumina/{species}_R1.fastq",
    r2 = "{species}/illumina/{species}_R2.fastq",
    long_reads = "{species}/nanopore/{species}.fastq",

    canu = "{species}/drafts/canu_{species}_contigs.fasta",
    flye = "{species}/drafts/flye_{species}_assembly.fasta",
    wengan = "{species}/drafts/wengan_{species}_assembly.fasta"
  params:
    root="{species}",
    busco_lineage=get_lineage
  output:
    "{species}/polished/best/polish.fasta"
  threads: 30
  script:
    "../scripts/polishing_pipeline.py"
