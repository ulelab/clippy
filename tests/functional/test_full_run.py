import clip
import os
import pybedtools
import pandas as pd
import cProfile


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


def test_single_gene_get_peaks_profiling(rootdir):
    xlinks = pybedtools.BedTool(
        os.path.join(rootdir, "tests", "data", "crosslinkcounts.bed")
    )
    annot = pd.read_table(
        os.path.join(rootdir, "tests", "data", "annot.gff"),
        header=None,
        names=[
            "chrom",
            "source",
            "feature_type",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ],
        comment="#",
    )
    annot_gene = annot[annot.feature_type == "gene"]
    annot_exons = annot[annot.feature_type == "exon"].copy(True)
    annot_exons["gene_id"] = [
        clip.get_gtf_attr_dict(attr_str)["gene_id"]
        for attr_str in annot_exons["attributes"]
    ]
    ang = pybedtools.BedTool.from_dataframe(annot_gene).sort()
    goverlaps = (
        xlinks.intersect(ang, s=True, wo=True)
        .to_dataframe(
            names=[
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "chrom2",
                "source",
                "feature",
                "gene_start",
                "gene_stop",
                "nothing",
                "strand2",
                "nothing2",
                "gene_name",
                "interval",
            ]
        )
        .drop(
            [
                "name",
                "chrom2",
                "nothing",
                "nothing2",
                "interval",
                "strand2",
                "source",
                "feature",
            ],
            axis=1,
        )
    )
    gene_annot_dict = {x: y for x, y in annot_gene.groupby("attributes")}
    arguments_list = [
        (
            pd.DataFrame(y),
            50,
            1,
            0.8,
            5,
            5,
            clip.get_exon_annot(x, annot_exons),
            None,
            gene_annot_dict[x],
        )
        for x, y in goverlaps.groupby("gene_name", as_index=False)
    ]
    pr = cProfile.Profile()
    pr.enable()
    output = [clip.single_gene_get_peaks(*args) for args in arguments_list[:1000]]
    pr.disable()
    pr.dump_stats(os.path.join(rootdir, "prof", "get_the_peaks.out"))
    assert len(output) == 1000


def test_getAllPeaks_yeast(rootdir, tmp_path):
    pr = cProfile.Profile()
    pr.enable()
    clip.getAllPeaks(
        pybedtools.BedTool(
            os.path.join(rootdir, "tests", "data", "crosslinkcounts.bed")
        ),
        os.path.join(rootdir, "tests", "data", "annot.gff"),
        50,
        1.0,
        0.8,
        5,
        5,
        1,
        1,
        os.path.join(tmp_path, "clippy_getAllPeaks_profile.bed"),
        False,
        None,
        0,
        0,
        os.path.join(rootdir, "tests", "data", "genome.fa.fai"),
        0,
    )
    pr.disable()
    pr.dump_stats(os.path.join(rootdir, "prof", "getAllPeaks.out"))


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
