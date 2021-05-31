import subprocess
from pathlib import Path

class BuscoResult:
  def __init__(self, *, contigs, busco_path, busco_score):
    self.contigs = contigs
    self.busco_score = busco_score
    self.busco_path = busco_path

class PolishPipeline:

  def __init__(self, *, root_dir, drafts, long_reads, short_reads, threads):
    self.root_dir = root_dir
    self.drafts = drafts
    self.long_reads = long_reads
    self.short_reads = short_reads
    self.MEDAKA_ROUNDS = 4
    self.PILON_ROUNDS = 4
    self.threads = threads

  def run(self):
    busco_results = [self.polish(draft) for draft in self.drafts]
    best_polish = max(busco_results, key=lambda busco_result: busco_result.busco_score)
    best_dir = f"{self.root_dir}/best"
    shell("mkdir -p {best_dir}")
    shell("cp {best_polish.contigs} {best_dir}")

  def polish(self, draft):
    polish_dir = f"{self.root_dir}/{draft.stem}"
    Path(polish_dir).mkdir(exist_ok=True)

    #pilon first
    pilon_polish = self.pilon_polish(polish_dir, draft)

    #then medaka
    medaka_polish = self.medaka_polish(polish_dir, pilon_polish.contigs)

    return medaka_polish


  def medaka_polish(self, polish_dir, draft):
    if not draft or not Path(draft).is_file():
      raise Exception(f"unknown draft type {draft}")

    # start with a new draft
    if draft not hasattr(draft, 'busco_score'):
      busco_result = BuscoResult(contigs=draft, busco_score=None, busco_path=None)
    # or a polished draft from pilon or other polishing tool
    else:
      busco_result = draft

    for i in range(self.MEDAKA_ROUNDS):
      out_dir = f"{polish_dir}/medaka/round_{i}"
      command = f"medaka_consensus -i {self.long_reads} -d {busco_result.contigs} -o {out_dir} -t {self.threads}"
      print(f"\n Running: {command} \n")
      subprocess.run(command.split(" "))

      polish = f"{out_dir}/consensus.fasta"
      if not Path(polish).is_file():
        raise Exception(f"medaka was unable to polish {busco_result.contigs} on round {i}")

      new_busco_result = run_busco(polish, f"{out_dir}/busco_out", "basidiomycota_odb10")

      #just ran busco for the best time
      if None in (busco_result.busco_score, busco_result.busco_path):
        busco_result = new_busco_result
      else:
        busco_result = determine_best_polish(new_busco_result, busco_result)

      if not busco_result:
        break

      return busco_result


  def pilon_polish(self, polish_dir, draft):
    if not draft or not Path(draft).is_file():
      raise Exception(f"unknown draft type {draft}")

    out_dir = f"{polish_dir}/pilon"
    Path(out_dir).mkdir(exist_ok=True)

    r1, r2 = self.short_reads

    # setting up inital busco results
    busco_result = BuscoResult(contigs=draft, busco_score=None, busco_path=None)
    for i in range(self.PILON_ROUNDS):
      sorted_aln = f"{out_dir}/sorted.bam"
      pilon_out = f"{out_dir}/round_{i}"
      print(f"\n ----------- {busco_result.contigs} -------------\n")

      shell("minimap2 -ax sr {busco_result.contigs} {r1} {r2} | samtools view -u | samtools sort -@ {self.threads} > {sorted_aln}")
      shell("samtools index {sorted_aln}")
      shell("pilon --genome {busco_result.contigs} --frags {sorted_aln} --outdir {pilon_out}")
      shell("rm {sorted_aln}")

      polish = f"{pilon_out}/pilon.fasta"
      print(f"\n ----------- {polish} -------------\n")

      if not Path(polish).is_file():
        raise Exception(f"pilon was unable to polish {busco_result.contigs} on round {i}")

      new_busco_result = run_busco(polish, f"{pilon_out}/busco_out", "basidiomycota_odb10")

      #just ran busco for the best time
      if None in (busco_result.busco_score, busco_result.busco_path):
        busco_result = new_busco_result
      else:
        busco_result = determine_best_polish(new_busco_result, busco_result)

      if not busco_result:
        break

    return busco_result


def determine_best_polish(new_busco_result, old_busco_result):
  if new_busco_result.busco_score > old_busco_result.busco_score:
    return new_busco_result
  else:
    return False


def get_busco_score(short_summary):
  """get busco Complete score from short_summary.txt"""

  with open(short_summary) as f:
    for line in f:
      line = line.strip()
      if line.startswith(("#","*")) or line == '':
        continue
      elif line.startswith("C:"):
        line = line.replace('%','').replace('[',',').replace(']','')
        return float(line.split(",")[0].split(':')[1])

# lineages/basidiomycota_odb10
def run_busco(contigs, outdir, lineage):
  busco_path = Path(outdir)

  #no slashes allowed in -o parameter so put stem as output
  cmd = f"busco -m genome -i {contigs} -o {busco_path.stem} -l {lineage} --cpu 25"
  print(f"\n ---------------------- {cmd} -------------------- \n")
  subprocess.run(cmd.split())
  shell("mv {busco_path.stem} {busco_path}")

  path = list(Path(outdir).glob("*/short_summary.txt"))
  if not path:
    raise Exception("cannot find short_summary.txt")

  short_summary = path[0]
  busco_score = get_busco_score(short_summary)
  return BuscoResult(
      contigs=contigs,
      busco_score=busco_score,
      busco_path=outdir
  )

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
        threads=30,
    )
    pipeline.run()
