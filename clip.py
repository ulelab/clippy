import numpy as np
import scipy.signal as sig
from scipy.ndimage.filters import uniform_filter1d
from scipy import stats
import pybedtools
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys

__version__ = "0.0.1"

def main():
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
    optional.add_argument('-mg',"--mingenecounts", type=int, default=5, nargs='?',
                        help='min counts per gene to look for peaks [DEFAULT 5]')
    optional.add_argument('-mb',"--minpeakcounts", type=int, default=5, nargs='?',
                        help='min counts per broad peak [DEFAULT 5]')
    optional.add_argument('-d',"--distance", type=int, default=10, nargs='?',
                        help='distance to merge crosslinks around single nt peaks [DEFAULT 10]')
    optional.add_argument('-g',"--mygene", type=str, nargs='?',
                        help='gene name, limits analysis to single gene')
    parser._action_groups.append(optional)
    args = parser.parse_args()
    print(args)
    outfile_name=args.outputprefix+"_rollmean" +str(args.windowsize,)+"_stdev"+str(args.adjust)+"_minGeneCount"+str(args.mingenecounts)+".bed"
    return(args.inputbed, args.annot, args.windowsize, args.adjust, args.mingenecounts, outfile_name, args.mygene, args.distance, args.minpeakcounts)

def getThePeaks(test, N, X, min_gene_count, counter):
    # Get the peaks for one gene
    # Now need to get an array of values
    chrom = test.chrom.iloc[0]
    genename = test.gene_name.iloc[0]
    strand = test.strand.iloc[0]
    start = int(test.gene_start.iloc[0])
    stop = int(test.gene_stop.iloc[0])
    all_vals = np.arange(start,stop,1)
    score_default = np.zeros(stop-start)
    default = pd.DataFrame({'start':all_vals, 'score':score_default})
    real = test.drop(['chrom','end','strand','gene_start','gene_stop','gene_name'],axis=1)
    mer = real.merge(default,how='outer') \
        .sort_values('score', ascending=False) \
        .drop_duplicates('start') \
        .sort_index() \
        .sort_values(by=['start'])

    scores = mer['score'].values
    if sum(scores) < min_gene_count:
        return(pd.DataFrame({'A' : []}), "", "")

    roll_mean_smoothed_scores = uniform_filter1d(scores.astype("float"), size=N)
    peaks=sig.find_peaks(roll_mean_smoothed_scores, height=np.mean(roll_mean_smoothed_scores), prominence=(np.std(roll_mean_smoothed_scores)*X))

    if peaks[0].size == 0:
        return(pd.DataFrame({'A' : []}), "", "")
    peaks_in_gene = []

    for i in range(0,len(peaks[0])):
        pk = peaks[0][i]
        final_peak = pd.DataFrame(data={"chrom":[chrom],"start":pk+start,"end":pk+start+1,"name":[genename],"score":["."],"strand":[strand]})
        peaks_in_gene.append(final_peak)

    peaks_in_gene = pd.concat(peaks_in_gene)
    peaks_in_gene = peaks_in_gene[["chrom","start","end","name","score","strand"]]

    if counter % 1000 == 0:
        print("Done for "+str(counter)+" genes")
    return(peaks_in_gene, roll_mean_smoothed_scores, peaks[0])

def getAllPeaks(counts_bed, annot, N, X, min_gene_count, outfile_name):
    pho92_iclip = pybedtools.BedTool(counts_bed)
    annot = pd.read_table(annot, header=None, names=["chrom","source","feature_type","start","end","score","strand","frame","attributes"], comment='#')
    annot_gene = annot[annot.feature_type=="gene"]
    ang = pybedtools.BedTool.from_dataframe(annot_gene).sort()
    goverlaps = pho92_iclip.intersect(ang, s=True, wo=True).to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand','chrom2','source','feature','gene_start', 'gene_stop','nothing','strand2','nothing2','gene_name','interval'])
    goverlaps.drop(['name','chrom2','nothing','nothing2','interval','strand2','source','feature'], axis=1, inplace=True)

    sep_genes = [pd.DataFrame(y) for x, y in goverlaps.groupby('gene_name', as_index=False)]

    all_peaks=[]
    counter =0
    for i in range(len(sep_genes)):
        df = sep_genes[i]
        counter += 1
        peaks_in_gene, rollingmean, plottingpeaks = getThePeaks(df, N, X, min_gene_count, counter)
        if peaks_in_gene.empty:
            continue
        else:
            all_peaks.append(peaks_in_gene)
    all_peaks = pd.concat(all_peaks)
    all_peaks.to_csv(outfile_name,sep="\t",header=False,index=False)
    return(pybedtools.BedTool.from_dataframe(all_peaks))
    print("Finished, written single nt peaks file.")

def getBroadPeaks(crosslinks, peaks, distance, min_peak_count, outfile_name): # crosslinks and peaks are both bedtools 
    merged_xlinks = crosslinks.merge(d=distance, c=[5,6], s=True, o=["sum","distinct"])
    final_peaks = merged_xlinks.intersect(peaks, s=True, u=True).filter(lambda x: float(x.score) >= min_peak_count)
    final_peaks.saveas(outfile_name)
    print("Finished, written broad peaks file.")


def getSingleGenePeaks(counts_bed, annot, N, X, min_gene_count, outfile_name, my_gene):
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

    peaks, roll_mean_smoothed_scores, plotting_peaks = getThePeaks(goverlaps, N, X, min_gene_count, 1)

    if peaks.empty:
        sys.exit("No peaks found in this gene with the current parameters")

    outfile_name=my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".bed"
    peaks.to_csv(outfile_name,sep="\t",header=False,index=False)

    # Make graph of gene
    plt.plot(roll_mean_smoothed_scores, '-bD', markevery=plotting_peaks.tolist())
    plt.ylabel('roll mean smoothed cDNAs')
    plt.axhline(y=np.mean(roll_mean_smoothed_scores),linewidth=4, color='r')
    plt.axhline(y=np.mean(roll_mean_smoothed_scores)+(np.std(roll_mean_smoothed_scores)*X),linewidth=1, color='g')
    plt.savefig(my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".png")
    print("Finished, written peak file and gene graph.")
    return(pybedtools.BedTool.from_dataframe(peaks))

if __name__ == "__main__":
    counts_bed, annot, N, X, min_gene_count, outfile_name, my_gene, distance, min_peak_count = main()
    counts_bed = pybedtools.BedTool(counts_bed)
    if not(my_gene is None):
        outfile_name=my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".bed"
        peaks = getSingleGenePeaks(counts_bed, annot, N, X, min_gene_count, outfile_name, my_gene)
        outfile_name=my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+"_broadPeaks.bed"
        getBroadPeaks(counts_bed, peaks, distance, min_peak_count, outfile_name)
    else:
        peaks = getAllPeaks(counts_bed, annot, N, X, min_gene_count, outfile_name)
        outfile_name=outfile_name.replace(".bed","_broadPeaks.bed")
        getBroadPeaks(counts_bed, peaks, distance, min_peak_count, outfile_name)

