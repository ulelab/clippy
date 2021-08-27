import pandas as pd
import pybedtools
import pybedtools.featurefuncs

genome_filepath = "tests/data/sk1_MvO_V1.fa.fai"

annot = pybedtools.BedTool("tests/data/annot.gff")

annot_bed = annot.filter(lambda row: row.fields[2] == "gene").saveas()

one_base_flanks = annot_bed.flank(
    g=genome_filepath,
    l=1,
    r=1,
    s=True
)

non_overlapping_flanks = one_base_flanks.intersect(annot_bed, v=True, s=True)

full_length_flanks = annot_bed.flank(
    g=genome_filepath,
    l=100,
    r=100,
    s=True
)

split_flanks = full_length_flanks.subtract(annot_bed, s=True)

def if_matching_return_a(interval):
    if interval.fields[8] == interval.fields[17]:
        return(
            pybedtools.cbedtools.create_interval_from_list(interval.fields[:9])
        )

valid_flanks = split_flanks \
    .intersect(non_overlapping_flanks, wo=True, s=True) \
    .each(if_matching_return_a) \
    .saveas()

assert(len(non_overlapping_flanks) == len(valid_flanks))

valid_flanks.each(pybedtools.featurefuncs.gff2bed).saveas("valid_flanks.bed")

