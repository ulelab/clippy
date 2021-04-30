import numpy as np
import scipy.signal as sig
from scipy.ndimage.filters import uniform_filter1d
from scipy import stats
import pybedtools
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import interaction

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
    optional.add_argument('-hc',"--height_cutoff", type=float, default=0.8, nargs='?',
                        help='proportion of prominence [DEFAULT 0.8]')
    optional.add_argument('-mg',"--mingenecounts", type=int, default=5, nargs='?',
                        help='min counts per gene to look for peaks [DEFAULT 5]')
    optional.add_argument('-mb',"--minpeakcounts", type=int, default=5, nargs='?',
                        help='min counts per broad peak [DEFAULT 5]')
    optional.add_argument('-g',"--mygene", type=str, nargs='?',
                        help='gene name, limits analysis to single gene')
    optional.add_argument('-int', "--interactive", action='store_true',
                        help='starts a Dash server to allow for interactive parameter tuning')
    parser._action_groups.append(optional)
    args = parser.parse_args()
    print(args)
    outfile_name=args.outputprefix+"_rollmean" +str(args.windowsize,)+"_stdev"+str(args.adjust)+"_minGeneCount"+str(args.mingenecounts)+".bed"
    return(args.inputbed, args.annot, args.windowsize, args.adjust, args.height_cutoff, args.mingenecounts,
        outfile_name, args.mygene, args.minpeakcounts, args.interactive)

def getThePeaks(test, N, X, rel_height, min_gene_count, counter):
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
    peaks=sig.find_peaks(roll_mean_smoothed_scores,
        height=np.mean(roll_mean_smoothed_scores),
        prominence=(np.std(roll_mean_smoothed_scores)*X),
        width=0.0,
        rel_height=rel_height)

    if peaks[0].size == 0:
        return(pd.DataFrame({'A' : []}), "", "")
    peaks_in_gene = []
    broad_peaks_in_gene = []
    for i in range(0,len(peaks[0])):
        pk = peaks[0][i]
        final_peak = pd.DataFrame(data={
            "chrom":[chrom],
            "start":pk+start,
            "end":pk+start+1,
            "name":[genename],
            "score":["."],
            "strand":[strand]
        })
        peaks_in_gene.append(final_peak)
        broad_peak = pd.DataFrame(data={
            "chrom":  [chrom],
            "start":  round(peaks[1]['left_ips'][i])+start,
            "end":    round(peaks[1]['right_ips'][i])+start+1,
            "name":   [genename],
            "score":  ["."],
            "strand": [strand]
        })
        broad_peaks_in_gene.append(broad_peak)

    peaks_in_gene = pd.concat(peaks_in_gene)
    peaks_in_gene = peaks_in_gene[["chrom","start","end","name","score","strand"]]

    broad_peaks_in_gene = pd.concat(broad_peaks_in_gene)
    broad_peaks_in_gene = broad_peaks_in_gene[["chrom","start","end","name","score","strand"]]

    if counter % 1000 == 0:
        print("Done for "+str(counter)+" genes")
    return(peaks_in_gene, broad_peaks_in_gene, roll_mean_smoothed_scores, peaks)

def getAllPeaks(counts_bed, annot, N, X, rel_height, min_gene_count, outfile_name):
    pho92_iclip = pybedtools.BedTool(counts_bed)
    annot = pd.read_table(annot, header=None, names=["chrom","source","feature_type","start","end","score","strand","frame","attributes"], comment='#')
    annot_gene = annot[annot.feature_type=="gene"]
    ang = pybedtools.BedTool.from_dataframe(annot_gene).sort()
    goverlaps = pho92_iclip.intersect(ang, s=True, wo=True).to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand','chrom2','source','feature','gene_start', 'gene_stop','nothing','strand2','nothing2','gene_name','interval'])
    goverlaps.drop(['name','chrom2','nothing','nothing2','interval','strand2','source','feature'], axis=1, inplace=True)

    sep_genes = [pd.DataFrame(y) for x, y in goverlaps.groupby('gene_name', as_index=False)]

    all_peaks=[]
    broad_peaks=[]
    counter =0
    for i in range(len(sep_genes)):
        df = sep_genes[i]
        counter += 1
        peaks_in_gene, broad_peaks_in_gene, rollingmean, peak_details = getThePeaks(df, N, X, rel_height, min_gene_count, counter)
        if peaks_in_gene.empty:
            continue
        else:
            all_peaks.append(peaks_in_gene)
            broad_peaks.append(broad_peaks_in_gene)
    all_peaks = pd.concat(all_peaks)
    all_peaks.to_csv(outfile_name,sep="\t",header=False,index=False)
    broad_peaks = pd.concat(broad_peaks)

    print("Finished, written single nt peaks file.")
    return(
        pybedtools.BedTool.from_dataframe(all_peaks),
        pybedtools.BedTool.from_dataframe(broad_peaks)
    )

def getBroadPeaks(crosslinks, broad_peaks, min_peak_count, outfile_name): # crosslinks and peaks are both bedtools 
    final_peaks = broad_peaks \
        .intersect(crosslinks, s=True, wo=True) \
        .sort() \
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

    peaks, broad_peaks, roll_mean_smoothed_scores, peak_details = getThePeaks(goverlaps, N, X, rel_height, min_gene_count, 1)

    if peaks.empty:
        sys.exit("No peaks found in this gene with the current parameters")

    outfile_name=my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".bed"
    peaks.to_csv(outfile_name,sep="\t",header=False,index=False)

    # Make graph of gene
    plt.plot(roll_mean_smoothed_scores, '-bD', markevery=peak_details[0].tolist())
    plt.ylabel('roll mean smoothed cDNAs')
    plt.axhline(y=np.mean(roll_mean_smoothed_scores),linewidth=4, color='r')
    plt.axhline(y=np.mean(roll_mean_smoothed_scores)+(np.std(roll_mean_smoothed_scores)*X),linewidth=1, color='g')
    plt.savefig(my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".png")
    print("Finished, written peak file and gene graph.")
    return(
        pybedtools.BedTool.from_dataframe(peaks),
        pybedtools.BedTool.from_dataframe(broad_peaks)
    )

if __name__ == "__main__":
    (counts_bed, annot, N, X, rel_height, min_gene_count, outfile_name, my_gene,
        min_peak_count, interactive) = main()
    counts_bed = pybedtools.BedTool(counts_bed)
    if interactive:
        app = interaction.DashApp(counts_bed, annot)
        app.run()
    else:
        if my_gene is None:
            peaks, broad_peaks = getAllPeaks(counts_bed, annot, N, X, rel_height, min_gene_count, outfile_name)
            outfile_name=outfile_name.replace(".bed","_broadPeaks.bed")
            getBroadPeaks(counts_bed, broad_peaks, min_peak_count, outfile_name)
        else:
            outfile_name=my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".bed"
            peaks, broad_peaks = getSingleGenePeaks(counts_bed, annot, N, X, rel_height, min_gene_count, outfile_name, my_gene)
            outfile_name=my_gene+"_rollmean" +str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+"_broadPeaks.bed"
            getBroadPeaks(counts_bed, broad_peaks, min_peak_count, outfile_name)
