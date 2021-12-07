import clip
import os


def test_alt_features_full_run(rootdir, monkeypatch, tmp_path):
    monkeypatch.setattr(
        "sys.argv",
        [
            "clippy",
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


def test_no_alt_features_cep295(rootdir, monkeypatch, tmp_path):
    monkeypatch.setattr(
        "sys.argv",
        [
            "clippy",
            "-i",
            os.path.join(rootdir, "tests", "data", "cep295.bed"),
            "-a",
            os.path.join(rootdir, "tests", "data", "gencode.v38.cep295.gtf"),
            "-o",
            os.path.join(tmp_path, "cep295"),
            "-g",
            os.path.join(rootdir, "tests", "data", "human_genome.fa.fai"),
        ],
    )
    clip.main()
    monkeypatch.setattr(
        "sys.argv",
        [
            "clippy",
            "-i",
            os.path.join(rootdir, "tests", "data", "cep295.bed"),
            "-a",
            os.path.join(rootdir, "tests", "data", "gencode.v38.cep295.gtf"),
            "-o",
            os.path.join(tmp_path, "cep295_scarna"),
            "-g",
            os.path.join(rootdir, "tests", "data", "human_genome.fa.fai"),
            "-alt",
            "scaRNA-gene_name-SCARNA[0-9]+",
        ],
    )
    clip.main()
    for file_name in os.listdir(tmp_path):
        test_lines = open(os.path.join(tmp_path, file_name)).readlines()
        ref_lines = open(
            os.path.join(rootdir, "tests", "data", "output", file_name)
        ).readlines()
        assert test_lines == ref_lines
