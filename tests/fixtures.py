"""Small, deterministic UniAnn inputs with two independent splice choices."""

from pathlib import Path


def write_fixture(directory: Path, splice_shift: int = 0) -> dict[str, Path]:
    """Write oriented scores and sequence; all positions in inputs are zero based.

    Each gene has a strong coding region and two acceptors three bases apart.
    Both introns preserve frame and exceed the model's 40-base minimum.
    The loci are separated so a 50-base flank cannot combine them.
    """
    directory.mkdir(parents=True, exist_ok=True)
    length = 1200
    sequence = list("C" * length)
    motifs = {
        "atg": ([60, 660], "ATG"),
        "stop": ([450, 1050], "TAA"),
        "gt": ([180 + splice_shift, 780 + splice_shift], "GT"),
        "ag": ([p + splice_shift for p in (238, 241, 838, 841)], "AG"),
    }
    paths = {name: directory / f"{name}.txt" for name in motifs}
    for name, (positions, motif) in motifs.items():
        for position in positions:
            sequence[position:position + len(motif)] = motif
        paths[name].write_text("".join(f"{p}\t2\n" for p in positions))
    paths["fasta"] = directory / "fixture.fa"
    paths["fasta"].write_text(">fixture\n" + "".join(sequence) + "\n")
    paths["emissions"] = directory / "emissions.txt"
    with paths["emissions"].open("w") as stream:
        for position in range(length):
            coding = any(start <= position < start + 393 for start in (60, 660))
            intronic = any(start + splice_shift <= position < start + splice_shift + 60
                           for start in (180, 780))
            values = [-8, -8, -8, -8, -8]
            values[0 if not coding else 1] = 2
            if intronic:
                values[4] = 3
            stream.write(str(position) + "\t" + "\t".join(map(str, values)) + "\n")
    return paths


def ordered_inputs(paths: dict[str, Path]) -> list[str]:
    return [str(paths[key]) for key in ("fasta", "emissions", "gt", "ag", "atg", "stop")]
