import clip
import os


def test_alt_features_full_run(rootdir, monkeypatch, tmp_path):
    monkeypatch.setattr(
        "sys.argv",
        [
            "clip.py",
            "-i",
            os.path.join(rootdir, "tests", "data", "crosslinkcounts.bed"),
            "-a",
            os.path.join(rootdir, "tests", "data", "annot.gff"),
            "-o",
            os.path.join(tmp_path, "clippy"),
            "-g",
            os.path.join(rootdir, "tests", "data", "genome.fa.fai"),
            "-t",
            str(1),
            "--alt_features",
            "scaRNA-gene_type-scaRNA,scaRNA-gene_name-SCARNA[0-9]+",
        ],
    )
    clip.main()
