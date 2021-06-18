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
            "-t",
            str(os.cpu_count()),
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
            "-t",
            str(1),
        ],
    )
    clip.main()
    for file_name in os.listdir(tmp_path):
        test_lines = open(os.path.join(tmp_path, file_name)).readlines()
        ref_lines = open(
            os.path.join(rootdir, "tests", "data", "output", file_name)
        ).readlines()
        assert test_lines == ref_lines


def test_getThePeaks_profiling(rootdir):
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
    arguments_list = [
        (pd.DataFrame(y), 50, 1, 0.8, 5, clip.get_exon_annot(x, annot_exons), None)
        for x, y in goverlaps.groupby("gene_name", as_index=False)
    ]
    pr = cProfile.Profile()
    pr.enable()
    output = [clip.getThePeaks(*args) for args in arguments_list[:1000]]
    pr.disable()
    pr.dump_stats(os.path.join(rootdir, "prof", "get_the_peaks.out"))
    assert len(output) == 1000
