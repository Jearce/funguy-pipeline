class Polish:
  def __init__(self, *, polish_path, busco_path, busco_score):
    self.polish_path = polish_path
    self.busco_score = busco_score
    self.busco_path = busco_path

class PolishPipeline:

  def __init__(self, *, root_dir, drafts, long_reads, short_reads):
    self.root_dir = root_dir
    self.drafts = drafts
    self.long_reads = long_reads
    self.short_reads = short_reads
    self.MEDAKA_ROUNDS = 4
    self.PILON_ROUNDS = 4

  def run(self):
    for index, draft in enumerate(self.drafts):
      polish_dir = f"{self.root_dir}/{draft.stem}"
      Path(polish_dir).mkdir(exist_ok=True)
      polished_assembly = self.pilon_polish(polish_dir, draft)

  def medaka_polish(self, polish_dir, draft):
    current_draft = draft
    for i in range(self.MEDAKA_ROUNDS):
      out_dir = f"{polish_dir}/medaka/round_{i}"
      command = f"medaka_consensus -i {self.long_reads} -d {current_draft} -o {out_dir} -t 10"
      print(f"\n Running: {command} \n")
      subprocess.run(command.split(" "))
      polish = Path(current_draft)
      if not polish.is_file():
        raise Exception(f"medaka was unable to polish {draft} on round {i}")

      #TODO: run busco and keep track of buso scores

      current_draft = f"{out_dir}/consensus.fasta"

  def pilon_polish(self, polish_dir, draft):
    out_dir = f"{polish_dir}/pilon"
    Path(out_dir).mkdir(exist_ok=True)

    r1, r2 = self.short_reads
    current_draft = draft
    for i in range(self.PILON_ROUNDS):
      sorted_aln = f"{out_dir}/sorted.bam"
      pilon_out = f"{out_dir}/round_{i}"

      shell("minimap2 -ax sr {current_draft} {r1} {r2} | samtools view -u | samtools sort -@ 10 > {sorted_aln}")
      shell("samtools index {sorted_aln}")
      shell("pilon --genome {current_draft} --frags {sorted_aln} --outdir {pilon_out}")
      shell("rm {sorted_aln}")

      polish = Path(current_draft)
      if not polish.is_file():
        raise Exception(f"pilon was unable to polish {draft} on round {i}")

      #TODO: run busco and keep track of buso scores
      busco_score = get_busco_score(current_draft)
      if busco_score > current_busco_score:
        best_polish = current_draft
      else:
        break


      current_draft = f"{pilon_out}/pilon.fasta"

# lineages/basidiomycota_odb10
def run_busco(contigs, outdir, lineage):
  command = f"busco -m genome -i {contigs} -o {outdir} -l {lineage} --cpu 20"

def get_busco_score(short_summary):
  with open(short_summary) as f:
    for line in f:
      line = line.strip()
      if line.startswith(("#","*")) or line == '':
        continue
      elif line.startswith("C:"):
        line = line.replace('%','').replace('[',',').replace(']','')
        return float(line.split(",")[0].split(':')[1])

rule reciprocal_polishing:
  input:
    drafts_dir = "{species}/drafts",
    r1 = "{species}/illumina/{species}_R1.fastq",
    r2 = "{species}/illumina/{species}_R2.fastq",
    long_reads = "{species}/nanopore/{species}.fastq"
  params:
    root="{species}"
  output:
    directory("{species}/polished")
  run:
    root = Path(input.drafts_dir)
    if not root.exists():
      raise Exception(f"${input.drafts_dir} does not exist")

    root_polish_dir =  f"{params.root}/polished"
    Path(root_polish_dir).mkdir(exist_ok=True)

    draft_assemblies = (file for file in root.glob("*.fasta") if file.is_file())
    pipeline = PolishPipeline(
        root_dir=root_polish_dir,
        drafts=draft_assemblies,
        long_reads=input.long_reads,
        short_reads=[input.r1, input.r2],
    )
    pipeline.run()
