import subprocess

from pathlib import Path

from snakemake.shell import shell


# so editor does not complain
snakemake = globals()["snakemake"]


class BuscoResult:
    def __init__(self, *, contigs, busco_path, busco_score):
        self.contigs = contigs
        self.busco_score = busco_score
        self.busco_path = busco_path


def initialize_busco_result(initial_draft):
    if not initial_draft or not Path(initial_draft).is_file():
        raise Exception(f"unknown draft type {initial_draft}")

    # start polished draft from other polishing tool
    if isinstance(initial_draft, BuscoResult):
        return initial_draft

    # or new draft
    else:
        return BuscoResult(
            contigs=initial_draft,
            busco_score=None,
            busco_path=None,
        )


def is_improved(new_busco_result, old_busco_result):
    if new_busco_result.busco_score > old_busco_result.busco_score:
        return True
    else:
        return False


def is_first(busco_result):
    if None in (busco_result.busco_score, busco_result.busco_path):
        return True
    else:
        return False


def get_busco_score(short_summary):
    """get busco Complete score from short_summary.txt"""

    with open(short_summary) as f:
        for line in f:
            line = line.strip()
            if line.startswith(("#", "*")) or line == "":
                continue
            elif line.startswith("C:"):
                line = line.replace("%", "").replace("[", ",").replace("]", "")
                return float(line.split(",")[0].split(":")[1])


# lineages/basidiomycota_odb10
def run_busco(contigs, outdir, lineage):
    busco_path = Path(outdir)

    # no slashes allowed in -o parameter so put stem as output
    cmd = f"busco -m genome -i {contigs} -o {busco_path.stem} -l {lineage} --cpu 25"
    subprocess.run(cmd.split())
    shell("mv {busco_path.stem} {busco_path}")

    path = list(Path(outdir).glob("*/short_summary.txt"))
    if not path:
        raise Exception("cannot find short_summary.txt")

    short_summary = path[0]
    busco_score = get_busco_score(short_summary)
    return BuscoResult(
        contigs=contigs, busco_score=busco_score, busco_path=outdir
    )


class PolishPipeline:
    def __init__(
        self, *, root_dir, drafts, long_reads, short_reads, threads, lineage
    ):
        self.root_dir = root_dir
        self.drafts = drafts
        self.long_reads = long_reads
        self.r1, self.r2 = short_reads
        self.MEDAKA_ROUNDS = 4
        self.PILON_ROUNDS = 4
        self.threads = threads
        self.busco_lineage = lineage

    def run(self):
        busco_results = [self.polish(draft) for draft in self.drafts]
        best_polish = max(
            busco_results, key=lambda busco_result: busco_result.busco_score
        )
        best_dir = f"{self.root_dir}/best"
        shell("mkdir -p {best_dir}")
        shell("cp {best_polish.contigs} {best_dir}/polish.fasta")

    def polish(self, draft):
        polish_dir = f"{self.root_dir}/{draft.stem}"
        Path(polish_dir).mkdir(exist_ok=True)

        # pilon first
        pilon_polish = self.pilon_polish(polish_dir, draft)

        # then medaka
        medaka_polish = self.medaka_polish(polish_dir, pilon_polish.contigs)

        return medaka_polish

    def medaka_polish(self, polish_dir, draft):
        busco_result = initialize_busco_result(draft)

        for i in range(self.MEDAKA_ROUNDS):
            out_dir = f"{polish_dir}/medaka/round_{i}"
            command = f"medaka_consensus -i {self.long_reads} -d {busco_result.contigs} -o {out_dir} -t {self.threads} -m r941_min_fast_g303"
            subprocess.run(command.split(" "))

            polish = f"{out_dir}/consensus.fasta"
            if not Path(polish).is_file():
                raise Exception(
                    f"medaka was unable to polish {busco_result.contigs} on round {i}"
                )

            new_busco_result = run_busco(
                polish, f"{out_dir}/busco_out", self.busco_lineage
            )

            # just ran first polish or the nth polish has improved assembly
            if is_first(busco_result) or is_improved(
                new_busco_result, busco_result
            ):
                busco_result = new_busco_result
            else:
                break

            return busco_result

    def pilon_polish(self, polish_dir, draft):
        busco_result = initialize_busco_result(draft)

        for i in range(self.PILON_ROUNDS):
            sorted_aln = "sorted.bam"
            pilon_out = f"{polish_dir}/pilon/round_{i}"
            shell(
                "minimap2 -ax sr {busco_result.contigs} {self.r1} {self.r2} | samtools view -u | samtools sort -@ {self.threads} > {sorted_aln}"
            )
            shell("samtools index {sorted_aln}")
            shell(
                "pilon --genome {busco_result.contigs} --frags {sorted_aln} --threads {self.threads} --outdir {pilon_out}"
            )
            shell("rm {sorted_aln}")
            polish = f"{pilon_out}/pilon.fasta"

            if not Path(polish).is_file():
                raise Exception(
                    f"pilon was unable to polish {busco_result.contigs} on round {i}"
                )

            new_busco_result = run_busco(
                polish, f"{pilon_out}/busco_out", self.busco_lineage
            )

            # just ran busco for the best time
            if is_first(busco_result) or is_improved(
                new_busco_result, busco_result
            ):
                busco_result = new_busco_result
            else:
                break

        return busco_result


def main():
    drafts_dir = f"{snakemake.params.root}/drafts"
    root = Path(drafts_dir)
    if not root.exists():
        raise Exception(f"{drafts_dir} does not exist")

    root_polish_dir = f"{snakemake.params.root}/polished"
    Path(root_polish_dir).mkdir(exist_ok=True)

    draft_assemblies = (
        file for file in root.glob("*.fasta") if file.is_file()
    )

    pipeline = PolishPipeline(
        root_dir=root_polish_dir,
        drafts=draft_assemblies,
        long_reads=snakemake.input.long_reads,
        short_reads=[snakemake.input.r1, snakemake.input.r2],
        threads=30,
        lineage=snakemake.params.busco_lineage,
    )
    pipeline.run()


main()
