import clip
import os
import pybedtools
import pandas as pd


def test_full_run_multi_thread(rootdir, monkeypatch, tmp_path):
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
            str(os.cpu_count()),
            "--no_exon_info",
        ],
    )
    clip.main()
    for file_name in os.listdir(tmp_path):
        test_lines = open(os.path.join(tmp_path, file_name)).readlines()
        ref_lines = open(
            os.path.join(rootdir, "tests", "data", "output", file_name)
        ).readlines()
        assert test_lines == ref_lines


def test_full_run_single_thread(rootdir, monkeypatch, tmp_path):
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
            "--no_exon_info",
        ],
    )
    clip.main()
    for file_name in os.listdir(tmp_path):
        test_lines = open(os.path.join(tmp_path, file_name)).readlines()
        ref_lines = open(
            os.path.join(rootdir, "tests", "data", "output", file_name)
        ).readlines()
        assert test_lines == ref_lines


def test_full_run_intergenic(rootdir, monkeypatch, tmp_path):
    monkeypatch.setattr(
        "sys.argv",
        [
            "clip.py",
            "-i",
            os.path.join(rootdir, "tests", "data", "crosslinkcounts.bed"),
            "-a",
            os.path.join(rootdir, "tests", "data", "annot.gff"),
            "-o",
            os.path.join(tmp_path, "intergenic_50"),
            "-g",
            os.path.join(rootdir, "tests", "data", "genome.fa.fai"),
            "-t",
            str(os.cpu_count()),
            "-inter",
            "50",
        ],
    )
    clip.main()
    for file_name in os.listdir(tmp_path):
        test_lines = open(os.path.join(tmp_path, file_name)).readlines()
        ref_lines = open(
            os.path.join(rootdir, "tests", "data", "output", file_name)
        ).readlines()
        assert test_lines == ref_lines
