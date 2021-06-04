import numpy as np
import scipy.signal as sig
from scipy.ndimage.filters import uniform_filter1d
from scipy import stats
import pybedtools
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import clip.interaction
from multiprocessing import Pool
import re

__version__ = "0.0.1"

def main():
    (counts_bed, annot, N, X, rel_height, min_gene_count, outfile_name, my_gene,
        min_peak_count, threads, chunksize_factor, interactive, exclusion_search) = parse_arguments(sys.argv[1:])
    counts_bed = pybedtools.BedTool(counts_bed)
    if interactive:
        app = clip.interaction.DashApp(counts_bed, annot)
        app.run()
    else:
        if my_gene is None:
            peaks, broad_peaks = getAllPeaks(counts_bed, annot, N, X, rel_height, min_gene_count, threads, chunksize_factor, outfile_name, exclusion_search)
            outfile_name=outfile_name.replace(".bed","_broadPeaks.bed")
            getBroadPeaks(counts_bed, broad_peaks, min_peak_count, outfile_name)
        else:
            outfile_name=my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".bed"
            peaks, broad_peaks = getSingleGenePeaks(counts_bed, annot, N, X, rel_height, min_gene_count, outfile_name, my_gene)
            outfile_name=my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+"_broadPeaks.bed"
            getBroadPeaks(counts_bed, broad_peaks, min_peak_count, outfile_name)

def parse_arguments(input_arguments):
    parser = argparse.ArgumentParser(description='Call CLIP peaks.')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i',"--inputbed", type=str, required=True,
                        help='bed file containing crosslink counts at each position')
    required.add_argument('-o',"--outputprefix", type=str, required=True,
                        help='prefix for output files')
    required.add_argument('-a',"--annot", type=str, required=True,
                        help='gtf annotation file')
    optional.add_argument('-n',"--windowsize", type=int, default=50, nargs='?',
                        help='rolling mean window size [DEFAULT 50]')
    optional.add_argument('-x',"--adjust", type=float, default=1, nargs='?',
                        help='adjustment for prominence [DEFAULT 1]')
    optional.add_argument('-hc',"--height_cutoff", type=float, default=0.8, nargs='?',
                        help='proportion of prominence [DEFAULT 0.8]')
    optional.add_argument('-mg',"--mingenecounts", type=int, default=5, nargs='?',
                        help='min counts per gene to look for peaks [DEFAULT 5]')
    optional.add_argument('-mb',"--minpeakcounts", type=int, default=5, nargs='?',
                        help='min counts per broad peak [DEFAULT 5]')
    optional.add_argument('-g',"--mygene", type=str, nargs='?',
                        help='gene name, limits analysis to single gene')
    optional.add_argument('-t', "--threads", type=int, default=1, nargs='?',
                        help='number of threads to use')
    optional.add_argument('-cf', "--chunksize_factor", type=int, default=16, nargs='?',
                        help='A factor used to control the number of jobs given to a thread at a time. A larger number reduces the number of jobs per chunk. Only increase if you experience crashes [DEFAULT 4]')
    optional.add_argument('-int', "--interactive", action='store_true',
                        help='starts a Dash server to allow for interactive parameter tuning')
    optional.add_argument('-ex',"--exclusion_search", type=str, nargs='?',
                        help='A list of GTF features to set different height thresholds on in the comma-separated format <gtf_key>-<search_pattern>-<threshold>')
    parser._action_groups.append(optional)
    args = parser.parse_args(input_arguments)
    print(args)
    outfile_name = ''.join([
        args.outputprefix,
        "_rollmean", str(args.windowsize),
        "_stdev", str(args.adjust),
        "_minGeneCount", str(args.mingenecounts),
        ".bed"])
    return(args.inputbed, args.annot, args.windowsize, args.adjust, args.height_cutoff,
        args.mingenecounts, outfile_name, args.mygene, args.minpeakcounts, args.threads,
        args.chunksize_factor, args.interactive, args.exclusion_search)

def getThePeaks(test, N, X, rel_height, min_gene_count):
    # Get the peaks for one gene
    # Now need to get an array of values
    chrom, xlink_start, xlink_end, score, strand, start, stop, gene_name = test.iloc[0]
    # BEDTools recognises GTF files for the intersection, but we have to take 1 away here
    start = int(start-1)
    stop = int(stop)
    xlink_coverage = {pos: 0 for pos in range(start, stop)}
    start_list = list(test.start)
    score_list = list(test.score)
    for idx in range(len(start_list)):
        xlink_coverage[start_list[idx]] += score_list[idx]
    scores = np.array(list(xlink_coverage.values()))

    if sum(scores) < min_gene_count:
        return(None, None, None, None)

    roll_mean_smoothed_scores = uniform_filter1d(scores.astype("float"), size=N)
    peaks=sig.find_peaks(roll_mean_smoothed_scores,
        height=np.mean(roll_mean_smoothed_scores),
        prominence=(np.std(roll_mean_smoothed_scores)*X),
        width=0.0,
        rel_height=rel_height)

    peak_num = len(peaks[0])
    if peak_num == 0:
        return(None, None, None, None)

    peaks_in_gene = np.array([
        (
            chrom,
            peaks[0][i]+start,
            peaks[0][i]+start+1,
            gene_name,
            ".",
            strand
        ) for i in range(peak_num)
    ])

    broad_peaks_in_gene = np.array([
        (
            chrom,
            round(peaks[1]['left_ips'][i])+start,
            round(peaks[1]['right_ips'][i])+start+1,
            gene_name,
            ".",
            strand
        ) for i in range(peak_num)
    ])

    return(peaks_in_gene, broad_peaks_in_gene, roll_mean_smoothed_scores, peaks)

def calc_chunksize(n_workers, len_iterable, factor):
    """Calculate chunksize argument for Pool-methods.
    #https://stackoverflow.com/questions/53751050/python-multiprocessing-understanding-logic-behind-chunksize
    """
    chunksize, extra = divmod(len_iterable, n_workers * factor)
    if extra:
        chunksize += 1
    return(chunksize)

def get_the_peaks_single_arg(input_tuple):
    return(getThePeaks(*input_tuple))

def get_gtf_attr_dict(attr_str):
    p = re.compile('(\S*).+"(\S*)"')
    return({
        m.group(1): m.group(2)
        for m in [
            p.search(attr.strip())
            for attr in attr_str.split(";")
        ] if m
    })

def parse_exclusion_search(exclusion_search_str):
    output = []
    searches = exclusion_search_str.split(",")
    for search in searches:
        split_search = search.split("-")
        output.append({
            "key": split_search[0],
            "regex": "-".join(split_search[1:(len(split_search)-1)]),
            "mean": float(split_search[-1])
        })
    return(output)

def test_exclusion(attr_str, exclusion_list):
    gtf_dict = get_gtf_attr_dict(attr_str)
    for exclusion in exclusion_list:
        if exclusion["key"] in gtf_dict:
            if re.search(exclusion["regex"], gtf_dict[exclusion["key"]]):
                return(True)
    return(False)

def getAllPeaks(counts_bed, annot, N, X, rel_height, min_gene_count, threads, chunksize_factor, outfile_name, exclusion_search):
    if threads > 1:
        pool = Pool(threads)
    pho92_iclip = pybedtools.BedTool(counts_bed)
    annot = pd.read_table(annot, header=None, names=["chrom","source","feature_type","start","end","score","strand","frame","attributes"], comment='#')
    annot_gene = annot[annot.feature_type=="gene"]

    if exclusion_search:
        exclusions = parse_exclusion_search(exclusion_search)
        annot_bool = [test_exclusion(attr_str, exclusions)
            for attr_str in annot_gene["attributes"]]

    ang = pybedtools.BedTool.from_dataframe(annot_gene).sort()
    goverlaps = pho92_iclip.intersect(ang, s=True, wo=True).to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand','chrom2','source','feature','gene_start', 'gene_stop','nothing','strand2','nothing2','gene_name','interval'])
    goverlaps.drop(['name','chrom2','nothing','nothing2','interval','strand2','source','feature'], axis=1, inplace=True)

    if threads > 1:
        arguments_list = [
            (pd.DataFrame(y), N, X, rel_height, min_gene_count)
            for x, y in goverlaps.groupby('gene_name', as_index=False)
        ]
        chunk_size = calc_chunksize(threads, len(arguments_list), chunksize_factor)
        output_list = pool.imap(get_the_peaks_single_arg, arguments_list, chunk_size)
    else:
        output_list = [getThePeaks(pd.DataFrame(y), N, X, rel_height, min_gene_count)
            for x, y in goverlaps.groupby('gene_name', as_index=False)]

    all_peaks=[]
    broad_peaks=[]
    for output in output_list:
        peaks_in_gene, broad_peaks_in_gene, rollingmean, peak_details = output
        if isinstance(peaks_in_gene, np.ndarray):
            all_peaks.append(peaks_in_gene)
            broad_peaks.append(broad_peaks_in_gene)

    all_peaks = pd.DataFrame(np.concatenate(all_peaks))
    all_peaks.to_csv(outfile_name,sep="\t",header=False,index=False)
    broad_peaks = pd.DataFrame(np.concatenate(broad_peaks))

    all_peaks_bed = pybedtools.BedTool.from_dataframe(all_peaks)\
        .sort()
    all_peaks_bed.saveas(outfile_name)

    print("Finished, written single nt peaks file.")
    return(
        all_peaks_bed,
        pybedtools.BedTool.from_dataframe(broad_peaks)
    )

def getBroadPeaks(crosslinks, broad_peaks, min_peak_count, outfile_name): # crosslinks and peaks are both bedtools
    # First, merge all broadpeaks
    final_peaks = broad_peaks \
        .sort() \
        .merge(s=True, c=[5,6], o=["distinct"]*2) \
    # Then, intersect with the crosslinks,
    # merge back down (to sum up the crosslink counts), and filter
    final_peaks = final_peaks \
        .intersect(crosslinks, s=True, wo=True) \
        .merge(s=True, c=[11, 6], o=["sum","distinct"]) \
        .filter(lambda x: float(x.score) >= min_peak_count)
    final_peaks.saveas(outfile_name)
    print("Finished, written broad peaks file.")

def getSingleGenePeaks(counts_bed, annot, N, X, rel_height, min_gene_count, outfile_name, my_gene):
    pho92_iclip = counts_bed
    annot = pd.read_table(annot, header=None, names=["chrom","source","feature_type","start","end","score","strand","frame","attributes"])
    annot_gene = annot[annot.feature_type=="gene"]
    # Search for the gene we want
    ang = annot_gene[annot_gene.attributes.str.contains(my_gene, case=False)]

    if len(ang) == 0:
        sys.exit("Couldn't find your gene in the provided annotation")
    elif len(ang) > 1:
        sys.exit("Found more than one gene containing that name - could you be more specific?")

    ang = pybedtools.BedTool.from_dataframe(ang)
 
    goverlaps = pho92_iclip.intersect(ang, s=True, wo=True).to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand','chrom2','source','feature','gene_start', 'gene_stop','nothing','strand2','nothing2','gene_name','interval'])
    goverlaps.drop(['name','chrom2','nothing','nothing2','interval','strand2','source','feature'], axis=1, inplace=True)

    peaks, broad_peaks, roll_mean_smoothed_scores, peak_details = getThePeaks(goverlaps, N, X, rel_height, min_gene_count)

    if not isinstance(peaks, np.ndarray):
        sys.exit("No peaks found in this gene with the current parameters")

    outfile_name=my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".bed"
    pd.DataFrame(peaks).to_csv(outfile_name,sep="\t",header=False,index=False)

    # Make graph of gene
    plt.plot(roll_mean_smoothed_scores, '-bD', markevery=peak_details[0].tolist())
    plt.ylabel('roll mean smoothed cDNAs')
    plt.axhline(y=np.mean(roll_mean_smoothed_scores),linewidth=4, color='r')
    plt.axhline(y=np.mean(roll_mean_smoothed_scores)+(np.std(roll_mean_smoothed_scores)*X),linewidth=1, color='g')
    plt.savefig(my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".png")
    print("Finished, written peak file and gene graph.")
    return(
        pybedtools.BedTool.from_dataframe(pd.DataFrame(peaks)),
        pybedtools.BedTool.from_dataframe(pd.DataFrame(broad_peaks))
    )

if __name__ == "__main__":
    main()
