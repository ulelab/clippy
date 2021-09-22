import numpy as np
from numpy.matrixlib.defmatrix import _from_string
import scipy.signal as sig
from scipy.ndimage.filters import uniform_filter1d
from scipy import stats
import pybedtools
import pybedtools.featurefuncs
import pandas as pd
import matplotlib
import tempfile

# matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import argparse
import sys
import clip.interaction
from multiprocessing import Pool
import re

__version__ = "1.1.0"


def main():
    (
        counts_bed,
        annot,
        N,
        X,
        rel_height,
        min_gene_count,
        outfile_name,
        my_gene,
        min_peak_count,
        threads,
        chunksize_factor,
        interactive,
        no_exon_info,
        alt_features,
        up_ext,
        down_ext,
        genome_file,
        intergenic_peak_threshold,
    ) = parse_arguments(sys.argv[1:])
    counts_bed = pybedtools.BedTool(counts_bed)
    if interactive:
        app = clip.interaction.DashApp(counts_bed, annot, genome_file)
        app.run()
    else:
        if my_gene is None:
            peaks, broad_peaks = getAllPeaks(
                counts_bed,
                annot,
                N,
                X,
                rel_height,
                min_gene_count,
                min_peak_count,
                threads,
                chunksize_factor,
                outfile_name,
                no_exon_info,
                alt_features,
                up_ext,
                down_ext,
                genome_file,
                intergenic_peak_threshold,
            )
            outfile_name = outfile_name.replace(".bed", "_broadPeaks.bed")
            getBroadPeaks(counts_bed, broad_peaks, min_peak_count, outfile_name)
        else:
            outfile_name = (
                my_gene
                + "_rollmean"
                + str(N)
                + "_stdev"
                + str(X)
                + "_minGeneCount"
                + str(min_gene_count)
                + ".bed"
            )
            peaks, broad_peaks = getSingleGenePeaks(
                counts_bed,
                annot,
                N,
                X,
                rel_height,
                min_gene_count,
                outfile_name,
                my_gene,
            )
            outfile_name = (
                my_gene
                + "_rollmean"
                + str(N)
                + "_stdev"
                + str(X)
                + "_minGeneCount"
                + str(min_gene_count)
                + "_broadPeaks.bed"
            )
            getBroadPeaks(counts_bed, broad_peaks, min_peak_count, outfile_name)


def parse_arguments(input_arguments):
    parser = argparse.ArgumentParser(description="Call CLIP peaks.")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i",
        "--inputbed",
        type=str,
        required=True,
        help="bed file containing crosslink counts at each position",
    )
    required.add_argument(
        "-o", "--outputprefix", type=str, required=True, help="prefix for output files",
    )
    required.add_argument(
        "-a", "--annot", type=str, required=True, help="gtf annotation file",
    )
    required.add_argument(
        "-g",
        "--genome_file",
        type=str,
        required=True,
        help=(
            "genome file containing chromosome lengths, used by BEDTools for "
            "genomic operations"
        ),
    )
    optional.add_argument(
        "-n",
        "--windowsize",
        type=int,
        default=50,
        nargs="?",
        help="rolling mean window size [DEFAULT 50]",
    )
    optional.add_argument(
        "-up",
        "--upstream_extension",
        type=int,
        default=0,
        nargs="?",
        help="upstream extension added to gene models [DEFAULT 0]",
    )
    optional.add_argument(
        "-down",
        "--downstream_extension",
        type=int,
        default=0,
        nargs="?",
        help="downstream extension added to gene models [DEFAULT 0]",
    )
    optional.add_argument(
        "-x",
        "--adjust",
        type=float,
        default=1,
        nargs="?",
        help="adjustment for prominence [DEFAULT 1]",
    )
    optional.add_argument(
        "-hc",
        "--height_cutoff",
        type=float,
        default=0.8,
        nargs="?",
        help="proportion of prominence [DEFAULT 0.8]",
    )
    optional.add_argument(
        "-mg",
        "--mingenecounts",
        type=int,
        default=5,
        nargs="?",
        help="min counts per gene to look for peaks [DEFAULT 5]",
    )
    optional.add_argument(
        "-mb",
        "--minpeakcounts",
        type=int,
        default=5,
        nargs="?",
        help="min counts per broad peak [DEFAULT 5]",
    )
    optional.add_argument(
        "-m",
        "--mygene",
        type=str,
        nargs="?",
        help="gene name, limits analysis to single gene",
    )
    optional.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        nargs="?",
        help="number of threads to use",
    )
    optional.add_argument(
        "-cf",
        "--chunksize_factor",
        type=int,
        default=16,
        nargs="?",
        help=(
            "A factor used to control the number of jobs given to a thread at "
            "a time. A larger number reduces the number of jobs per chunk. "
            "Only increase if you experience crashes [DEFAULT 4]"
        ),
    )
    optional.add_argument(
        "-int",
        "--interactive",
        action="store_true",
        help="starts a Dash server to allow for interactive parameter tuning",
    )
    optional.add_argument(
        "-nei",
        "--no_exon_info",
        action="store_true",
        help="Turn off individual exon and intron thresholds",
    )
    optional.add_argument(
        "-inter",
        "--intergenic_peak_threshold",
        type=int,
        default=0,
        nargs="?",
        help=(
            "Intergenic peaks are called by first creating intergenic regions "
            "and calling peaks on the regions as though they were genes. "
            "The regions are made by expanding intergenic crosslinks and "
            "merging the result. This parameter is the threshold number of "
            "crosslinks required to include a region. If set to zero (default), "
            "no intergenic peaks will be called. When using this mode, the "
            "intergenic regions used will be output as a GTF file. [DEFAULT 0]"
        ),
    )
    optional.add_argument(
        "-alt",
        "--alt_features",
        type=str,
        nargs="?",
        help=(
            "A list of alternative GTF features to set individual height "
            "thresholds on in the comma-separated format "
            "<alt_feature_name>-<gtf_key>-<search_pattern>"
        ),
    )
    parser._action_groups.append(optional)
    args = parser.parse_args(input_arguments)
    print(args)
    outfile_name = "".join(
        [
            args.outputprefix,
            "_rollmean",
            str(args.windowsize),
            "_stdev",
            str(args.adjust),
            "_minGeneCount",
            str(args.mingenecounts),
            ".bed",
        ]
    )
    return (
        args.inputbed,
        args.annot,
        args.windowsize,
        args.adjust,
        args.height_cutoff,
        args.mingenecounts,
        outfile_name,
        args.mygene,
        args.minpeakcounts,
        args.threads,
        args.chunksize_factor,
        args.interactive,
        args.no_exon_info,
        args.alt_features,
        args.upstream_extension,
        args.downstream_extension,
        args.genome_file,
        args.intergenic_peak_threshold,
    )


def single_gene_get_peaks(
    test,
    N,
    X,
    rel_height,
    min_gene_count,
    min_peak_count,
    annot_exon,
    alt_features,
    gene_flank_annot,
):
    # Get the peaks for one gene
    # Now need to get an array of values
    chrom, xlink_start, xlink_end, score, strand, start, stop, gene_name = test.iloc[0]
    # The gene flanks (if set) will be passed in, which we need to incorporate
    start = min(gene_flank_annot.start)
    stop = max(gene_flank_annot.end)
    # BEDTools recognises GTF files for the intersection, but we have to take 1 away here
    start = int(start - 1)
    stop = int(stop)
    xlink_coverage = {pos: 0 for pos in range(start, stop)}
    start_list = list(test.start)
    score_list = list(test.score)
    for idx in range(len(start_list)):
        xlink_coverage[start_list[idx]] += score_list[idx]
    scores = np.array(list(xlink_coverage.values()))

    if sum(scores) < min_gene_count:
        return (None, None, None, None, None, None)

    feature_names = set(["intron", "exon"])

    nucleotide_sets = [set(["intron"]) for i in range(stop - start)]

    if isinstance(annot_exon, pd.DataFrame):
        for index, row in annot_exon.iterrows():
            for pos in range(row["start"] - 1, row["end"]):
                zero_pos = pos - start
                if "intron" in nucleotide_sets[zero_pos]:
                    nucleotide_sets[zero_pos].remove("intron")
                nucleotide_sets[zero_pos].add("exon")

    if alt_features:
        for feature_name, alt_features_annot in alt_features.items():
            if len(alt_features_annot) > 0:
                threshold_overrides = (
                    pybedtools.BedTool.from_dataframe(alt_features_annot)
                    .intersect(
                        pybedtools.BedTool(
                            "\t".join(
                                [chrom, str(start), str(stop), "gene", ".", strand]
                            ),
                            from_string=True,
                        ),
                        s=True,
                    )
                    .to_dataframe()
                )
                for index, row in threshold_overrides.iterrows():
                    for pos in range(row["start"] - 1, row["end"]):
                        zero_pos = pos - start
                        for name in ["intron", "exon"]:
                            if name in nucleotide_sets[zero_pos]:
                                nucleotide_sets[zero_pos].remove(name)
                        nucleotide_sets[zero_pos].add(feature_name)
                        feature_names.add(feature_name)

    roll_mean_smoothed_scores = uniform_filter1d(scores.astype("float"), size=N)

    values_dict = {}
    for feature_name in feature_names:
        feature_values = roll_mean_smoothed_scores[
            [feature_name in i for i in nucleotide_sets]
        ]
        # Removes zeros from the vectors, but not sure if this is a good idea
        # Especially for calculating standard deviation
        # feature_values = feature_values[np.logical_not(feature_values == 0.0)]
        values_dict[feature_name] = feature_values

    mean_dict = {}
    std_dict = {}
    for feature_name in feature_names:
        if len(values_dict[feature_name]) > 0:
            mean_dict[feature_name] = np.mean(values_dict[feature_name])
            std_dict[feature_name] = np.std(values_dict[feature_name])
        else:
            mean_dict[feature_name] = 0.0
            std_dict[feature_name] = 0.0

    # In the case that the intron mean is greater than the exon mean, use the
    #  entire gene for the threshold
    if mean_dict["intron"] > mean_dict["exon"]:
        mean_dict["exon"] = np.mean(
            np.concatenate((values_dict["intron"], values_dict["exon"]))
        )
        std_dict["exon"] = np.std(
            np.concatenate((values_dict["intron"], values_dict["exon"]))
        )

    heights = np.array(
        [max([mean_dict[feature] for feature in pos]) for pos in nucleotide_sets]
    )

    prominences = (
        np.array(
            [max([std_dict[feature] for feature in pos]) for pos in nucleotide_sets]
        )
        * X
    )

    peaks = sig.find_peaks(
        roll_mean_smoothed_scores,
        height=heights,
        prominence=prominences,
        width=0.0,
        rel_height=rel_height,
    )

    peak_num = len(peaks[0])
    if peak_num == 0:
        return (None, None, None, None, None, None)

    peaks_in_gene = np.array(
        [
            (
                chrom,
                peaks[0][i] + start,
                peaks[0][i] + start + 1,
                gene_name,
                ".",
                strand,
            )
            for i in range(peak_num)
        ]
    )

    merged_broad_peaks = []

    for j in range(peak_num):
        new_peak = {
            "start": round(peaks[1]["left_ips"][j]) + start,
            "end": round(peaks[1]["right_ips"][j]) + start + 1,
        }
        added = False
        for peak in merged_broad_peaks:
            if (
                peak["start"] <= new_peak["start"] and new_peak["start"] <= peak["end"]
            ) or (peak["start"] <= new_peak["end"] and new_peak["end"] <= peak["end"]):
                peak["start"] = min(peak["start"], new_peak["start"])
                peak["end"] = max(peak["end"], new_peak["end"])
                added = True
        if not added:
            merged_broad_peaks.append(new_peak)

    broad_peaks_in_gene = np.array(
        [
            (
                chrom,
                peak["start"],
                peak["end"],
                gene_name,
                scores[peak["start"] - start : peak["end"] - start].sum(),
                strand,
            )
            for peak in merged_broad_peaks
        ]
    )

    filt_broad_peaks_in_gene = broad_peaks_in_gene[
        broad_peaks_in_gene[:, 4].astype("float") >= min_peak_count
    ]

    return (
        peaks_in_gene,
        filt_broad_peaks_in_gene,
        roll_mean_smoothed_scores,
        peaks,
        heights,
        prominences,
    )


def calc_chunksize(n_workers, len_iterable, factor):
    """Calculate chunksize argument for Pool-methods.
    https://stackoverflow.com/questions/53751050/python-multiprocessing-understanding-logic-behind-chunksize
    """
    chunksize, extra = divmod(len_iterable, n_workers * factor)
    if extra:
        chunksize += 1
    return chunksize


def get_the_peaks_single_arg(input_tuple):
    return single_gene_get_peaks(*input_tuple)


def get_gtf_attr_dict(attr_str):
    p = re.compile(r'(\S*).+"(\S*)"')
    return {
        m.group(1): m.group(2)
        for m in [p.search(attr.strip()) for attr in attr_str.split(";")]
        if m
    }


def parse_alt_features(alt_features_search_str):
    output = {}
    searches = alt_features_search_str.split(",")
    for search in searches:
        split_search = search.split("-")
        if len(split_search) >= 3:
            output.setdefault(split_search[0], [])
            output[split_search[0]].append(
                {"key": split_search[1], "regex": "-".join(split_search[2:]),}
            )
    return output


def test_alt_features(attr_str, alt_features_dict):
    return_set = set([])
    gtf_dict = get_gtf_attr_dict(attr_str)
    for alt_feature_name, alt_features_list in alt_features_dict.items():
        for alt_feature in alt_features_list:
            if alt_feature["key"] in gtf_dict:
                if re.search(alt_feature["regex"], gtf_dict[alt_feature["key"]]):
                    return_set.add(alt_feature_name)
    return return_set


def return_alt_features(alt_feature_str, annot_gene):
    alt_features = parse_alt_features(alt_feature_str)
    alt_feature_matches = [
        test_alt_features(attr_str, alt_features)
        for attr_str in annot_gene["attributes"]
    ]
    annot_alt_features = {
        feature_name: pd.DataFrame(
            annot_gene.iloc[[feature_name in i for i in alt_feature_matches]], copy=True
        ).reset_index(drop=True)
        for feature_name in alt_features.keys()
    }
    return annot_alt_features


def get_exon_annot(gene_attrs, annot_exons):
    if isinstance(annot_exons, dict):
        attr_dict = get_gtf_attr_dict(gene_attrs)
        if "gene_id" in attr_dict:
            if attr_dict["gene_id"] in annot_exons:
                return annot_exons[attr_dict["gene_id"]]
    return None

def get_overlapping_feature_bed(input_annot_bed, genome_file):
    # Find regions of the genome which are overlapped by multiple features on
    # the same strand and remove those crosslinks.
    overlapping_feature_dfs = []
    for strand in ["+", "-"]:
        genomecov = input_annot_bed.genome_coverage(
            bg=True, strand=strand, g=genome_file
        )
        if genomecov.count() > 0:
            genomecov = genomecov.to_dataframe()
            genomecov = genomecov[genomecov.name > 1].copy()
            genomecov["score"] = genomecov.name
            genomecov["strand"] = strand
            overlapping_feature_dfs.append(genomecov)
    if len(overlapping_feature_dfs) == 0:
        return None
    return pybedtools.BedTool.from_dataframe(pd.concat(overlapping_feature_dfs)).sort()


def get_flanking_regions(annot_gene_bed, genome_filepath, up_ext, down_ext):
    # Get bases either end of each gene
    one_base_flanks = annot_gene_bed.flank(g=genome_filepath, l=1, r=1, s=True)
    # Get the bases which don't overlap any other gene
    non_overlapping_one_base_flanks = one_base_flanks.intersect(
        annot_gene_bed, v=True, s=True
    )
    # Get the full length flanks
    full_length_flanks = annot_gene_bed.flank(
        g=genome_filepath, l=up_ext, r=down_ext, s=True
    )
    # Remove any parts of the full length flanks which overlap genes
    split_flanks = full_length_flanks.subtract(annot_gene_bed, s=True)
    # A function which takes a row of BEDTools intersect output (when run with
    # the -wo flag) and returns the first entry (from file a) if the attribute
    # fields match
    def if_matching_return_a(interval):
        if interval.fields[8] == interval.fields[17]:
            return pybedtools.cbedtools.create_interval_from_list(interval.fields[:9])

    #  Converts the gene field to flank, so that the resulting flanks can be
    # identified
    def convert_type_to_flank(interval):
        assert interval.file_type == "gff"
        new_fields = interval.fields[:]
        new_fields[2] = "flank"
        return pybedtools.cbedtools.create_interval_from_list(new_fields)

    # Intersect the split flanks with the one base flanks which don't overlap
    # genes and filter for matching attribute fields
    valid_flanks = (
        split_flanks.intersect(non_overlapping_one_base_flanks, wo=True, s=True)
        .each(if_matching_return_a)
        .each(convert_type_to_flank)
        .saveas()
    )
    # Assertion that the one valid flank is returned for every non-overlapping
    # one-base flank
    assert len(non_overlapping_one_base_flanks) == len(valid_flanks)
    return valid_flanks


def getAllPeaks(
    counts_bed,
    annot_filepath,
    N,
    X,
    rel_height,
    min_gene_count,
    min_peak_count,
    threads,
    chunksize_factor,
    outfile_name,
    no_exon_info,
    alt_features,
    up_ext,
    down_ext,
    genome_file,
    intergenic_peak_threshold,
):
    if threads > 1:
        pool = Pool(threads)
    clip_bed = pybedtools.BedTool(counts_bed)
    annot = pd.read_table(
        annot_filepath,
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
    annot_gene = annot[annot.feature_type == "gene"].copy(True)

    # Set up exon information
    annot_exons = None
    if not no_exon_info:
        annot_exons = annot[annot.feature_type == "exon"].copy(True)
        annot_exons["gene_id"] = [
            get_gtf_attr_dict(attr_str)["gene_id"]
            for attr_str in annot_exons["attributes"]
        ]
        annot_exons = {
            x: y for x, y in annot_exons.groupby("gene_id", as_index=False)
        }


    #  Setup alt_features
    annot_alt_features = None
    if alt_features:
        annot_alt_features = return_alt_features(alt_features, annot_gene)

    annot_bed = pybedtools.BedTool.from_dataframe(annot_gene).sort()

    flank_bed = get_flanking_regions(
        annot_bed, genome_file, max(N, up_ext), max(N, down_ext)
    )

    annot_bed_with_flanks = annot_bed.cat(flank_bed, postmerge=False).sort()

    # if intergenic_peak_threshold==0 then no intergenic peaks will be called
    if intergenic_peak_threshold > 0:

        def add_intergenic_gff_attribute_field(interval):
            new_fields = interval.fields[:]
            new_fields[2] = "intergenic_region"
            new_fields[8] = 'name="intergenic_region_{}_{}_{}_{}";'.format(
                new_fields[0], new_fields[3], new_fields[4], new_fields[6],
            )
            return pybedtools.cbedtools.create_interval_from_list(new_fields)

        # Get the xlinks which don't overlap genes or flanks, expand them,
        # subtract the genes and flanks (in case you expanded into them)
        # merge the result, filter for regions with enough crosslinks, and
        # convert to gff
        intergenic_regions = (
            clip_bed.intersect(annot_bed_with_flanks, s=True, v=True)
            .slop(g=genome_file, b=N)
            .subtract(annot_bed_with_flanks, s=True)
            .sort()
            .merge(s=True, c=[5, 6], o=["sum", "distinct"])
            .filter(lambda row: int(row.score) >= intergenic_peak_threshold)
            .each(pybedtools.featurefuncs.bed2gff)
            .each(add_intergenic_gff_attribute_field)
            .saveas(outfile_name.replace(".bed", "_intergenic_regions.gtf"))
        )
        annot_bed_with_flanks = annot_bed_with_flanks.cat(
            intergenic_regions, postmerge=False
        ).sort()

    # Split crosslinks based on overlap with features
    goverlaps = clip_bed.intersect(annot_bed_with_flanks, s=True, wo=True).to_dataframe(
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
    goverlaps.drop(
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
        inplace=True,
    )

    gene_flank_dict = {
        x: y
        for x, y in annot_bed_with_flanks.to_dataframe().groupby(
            "attributes", as_index=False
        )
    }

    if threads > 1:
        arguments_list = [
            (
                pd.DataFrame(y),
                N,
                X,
                rel_height,
                min_gene_count,
                min_peak_count,
                get_exon_annot(x, annot_exons),
                annot_alt_features,
                gene_flank_dict[x],
            )
            for x, y in goverlaps.groupby("gene_name", as_index=False)
        ]
        chunk_size = calc_chunksize(threads, len(arguments_list), chunksize_factor)
        output_list = pool.imap(get_the_peaks_single_arg, arguments_list, chunk_size)
    else:
        output_list = [
            single_gene_get_peaks(
                pd.DataFrame(y),
                N,
                X,
                rel_height,
                min_gene_count,
                min_peak_count,
                get_exon_annot(x, annot_exons),
                annot_alt_features,
                gene_flank_dict[x],
            )
            for x, y in goverlaps.groupby("gene_name", as_index=False)
        ]

    all_peaks = []
    broad_peaks = []
    for output in output_list:
        (
            peaks_in_gene,
            broad_peaks_in_gene,
            rollingmean,
            peak_details,
            heights,
            prominences,
        ) = output
        if isinstance(peaks_in_gene, np.ndarray):
            all_peaks.append(peaks_in_gene)
            broad_peaks.append(broad_peaks_in_gene)

    all_peaks = pd.DataFrame(np.concatenate(all_peaks))
    all_peaks.to_csv(outfile_name, sep="\t", header=False, index=False)
    broad_peaks = pd.DataFrame(np.concatenate(broad_peaks))

    overlapping_features = get_overlapping_feature_bed(
        annot_bed_with_flanks, genome_file
    )
    filtered_broad_peaks = None
    if overlapping_features is None:
        filtered_broad_peaks = pybedtools.BedTool.from_dataframe(broad_peaks)
    else:
        filtered_broad_peaks = pybedtools.BedTool.from_dataframe(broad_peaks).intersect(
            overlapping_features, v=True, s=True
        )

    all_peaks_bed = pybedtools.BedTool.from_dataframe(all_peaks).sort()
    all_peaks_bed.saveas(outfile_name)

    print("Finished, written single nt peaks file.")
    return (all_peaks_bed, filtered_broad_peaks)


def getBroadPeaks(
    crosslinks, broad_peaks, min_peak_count, outfile_name
):  # crosslinks and peaks are both bedtools
    # First, merge all broadpeaks
    final_peaks = broad_peaks.sort().merge(
        s=True, c=[5, 6], o=["distinct"] * 2
    )  # Then, intersect with the crosslinks,
    # merge back down (to sum up the crosslink counts), and filter
    final_peaks = (
        final_peaks.intersect(crosslinks, s=True, wo=True)
        .merge(s=True, c=[11, 6], o=["sum", "distinct"])
        .filter(lambda x: float(x.score) >= min_peak_count)
    )
    final_peaks.saveas(outfile_name)
    print("Finished, written broad peaks file.")


def getSingleGenePeaks(
    counts_bed,
    annot,
    N,
    X,
    rel_height,
    min_gene_count,
    outfile_name,
    my_gene,
    min_peak_count,
):
    pho92_iclip = counts_bed
    annot = pd.read_table(
        annot,
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
    )
    annot_gene = annot[annot.feature_type == "gene"]
    # Search for the gene we want
    ang = annot_gene[annot_gene.attributes.str.contains(my_gene, case=False)]

    if len(ang) == 0:
        sys.exit("Couldn't find your gene in the provided annotation")
    elif len(ang) > 1:
        sys.exit(
            "Found more than one gene containing that name "
            "- could you be more specific?"
        )

    ang = pybedtools.BedTool.from_dataframe(ang)

    goverlaps = pho92_iclip.intersect(ang, s=True, wo=True).to_dataframe(
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
    goverlaps.drop(
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
        inplace=True,
    )

    peaks, broad_peaks, roll_mean_smoothed_scores, peak_details = single_gene_get_peaks(
        goverlaps, N, X, rel_height, min_gene_count, min_peak_count
    )

    if not isinstance(peaks, np.ndarray):
        sys.exit("No peaks found in this gene with the current parameters")

    outfile_name = (
        my_gene
        + "_rollmean"
        + str(N)
        + "_stdev"
        + str(X)
        + "_minGeneCount"
        + str(min_gene_count)
        + ".bed"
    )
    pd.DataFrame(peaks).to_csv(outfile_name, sep="\t", header=False, index=False)

    # Make graph of gene
    plt.plot(roll_mean_smoothed_scores, "-bD", markevery=peak_details[0].tolist())
    plt.ylabel("roll mean smoothed cDNAs")
    plt.axhline(y=np.mean(roll_mean_smoothed_scores), linewidth=4, color="r")
    plt.axhline(
        y=np.mean(roll_mean_smoothed_scores) + (np.std(roll_mean_smoothed_scores) * X),
        linewidth=1,
        color="g",
    )
    plt.savefig(
        my_gene
        + "_rollmean"
        + str(N)
        + "_stdev"
        + str(X)
        + "_minGeneCount"
        + str(min_gene_count)
        + ".png"
    )
    print("Finished, written peak file and gene graph.")
    return (
        pybedtools.BedTool.from_dataframe(pd.DataFrame(peaks)),
        pybedtools.BedTool.from_dataframe(pd.DataFrame(broad_peaks)),
    )


if __name__ == "__main__":
    main()
