import numpy as np
import scipy.signal as sig
from scipy.ndimage.filters import uniform_filter1d
from scipy import stats
import pybedtools
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Call CLIP peaks.')
    parser.add_argument('-i',"--inputbed", type=str,
                        help='bed file containing crosslink counts at each position')
    parser.add_argument('-o',"--outputprefix", type=str,
                        help='prefix for output files')
    parser.add_argument('-a',"--annot", type=str,
                        help='gtf annotation file')
    parser.add_argument('-n',"--windowsize", type=int,
                        help='rolling mean window size')
    parser.add_argument('-x',"--adjust", type=int,
                        help='adjustment for prominence')
    parser.add_argument('-m',"--mingenecounts", type=int,
                        help='min counts per gene to look for peaks')
    parser.add_argument('-g',"--mygene", type=str,
                        help='gene name, limits analysis to single gene')
    args = parser.parse_args()
    outfile_name=args.outputprefix+str(args.windowsize)+"_stdev"+str(args.adjust)+"_minGeneCount"+str(args.mingenecounts)+".bed"
    return(args.inputbed, args.annot, args.windowsize, args.adjust, args.mingenecounts, outfile_name, args.mygene)

def getThePeaks(test, N, X, min_gene_count):
    # Get the peaks for one gene
    # Now need to get an array of values
    chrom = test.chrom.iloc[0]
    genename = test.gene_name.iloc[0]
    strand = test.strand.iloc[0]
    start = test.gene_start.iloc[0]
    stop = test.gene_stop.iloc[0]
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
        return pd.DataFrame({'A' : []})
    #print(sig.find_peaks(scores, height=np.median(scores)))
    roll_mean_smoothed_scores = uniform_filter1d(scores.astype("float"), size=N)
    peaks=sig.find_peaks(roll_mean_smoothed_scores, height=np.mean(roll_mean_smoothed_scores), prominence=(np.std(roll_mean_smoothed_scores)*X))
    if peaks[0].size == 0:
        return pd.DataFrame({'A' : []})
    peaks_in_gene = []

    for i in range(0,len(peaks[0])):
        pk = peaks[0][i]
        #peak_start = peaks[1]['left_bases'][i]
        #peak_end = peaks[1]['right_bases'][i]
        final_peak = pd.DataFrame(data={"chrom":[chrom],"start":pk+start,"end":pk+start+1,"name":[genename],"score":["."],"strand":[strand]})
        peaks_in_gene.append(final_peak)

    peaks_in_gene = pd.concat(peaks_in_gene)
    return(peaks_in_gene, roll_mean_smoothed_scores, peaks[0])

def getAllPeaks(counts_bed, annot, N, X, min_gene_count, outfile_name):
    pho92_iclip = pybedtools.BedTool(counts_bed)
    annot = pd.read_table(annot, header=None, names=["chrom","source","feature_type","start","end","score","strand","frame","attributes"])
    annot_gene = annot[annot.feature_type=="gene"]
    ang = pybedtools.BedTool.from_dataframe(annot_gene).sort()
    goverlaps = pho92_iclip.intersect(ang, s=True, wo=True).to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand','chrom2','source','feature','gene_start', 'gene_stop','nothing','strand2','nothing2','gene_name','interval'])
    goverlaps.drop(['name','chrom2','nothing','nothing2','interval','strand2','source','feature'], axis=1, inplace=True)

    sep_genes = [pd.DataFrame(y) for x, y in goverlaps.groupby('gene_name', as_index=False)]

    all_peaks=[]
    for df in sep_genes:
        peaks_in_gene = getThePeaks(df)[0]
        if peaks_in_gene.empty:
            continue
        else:
            all_peaks.append(peaks_in_gene)
    all_peaks = pd.concat(all_peaks)
    all_peaks.to_csv(outfile_name,sep="\t",header=False,index=False)
    print("Finished, written peaks file.")

def getSingleGenePeaks(counts_bed, annot, N, X, min_gene_count, outfile_name, my_gene):
    pho92_iclip = pybedtools.BedTool(counts_bed)
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

    peaks, roll_mean_smoothed_scores, plotting_peaks = getThePeaks(goverlaps, N, X, min_gene_count)

    if peaks.empty:
        sys.exit("No peaks found in this gene with the current parameters")

    outfile_name=my_gene+str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".bed"
    peaks.to_csv(outfile_name,sep="\t",header=False,index=False)

    # Make graph of gene
    print(roll_mean_smoothed_scores)
    print(peaks)
    plt.plot(roll_mean_smoothed_scores, '-bD', markevery=plotting_peaks.tolist())
    plt.ylabel('roll mean smoothed cDNAs')
    plt.axhline(y=np.mean(roll_mean_smoothed_scores),linewidth=4, color='r')
    plt.axhline(y=np.mean(roll_mean_smoothed_scores)+(np.std(roll_mean_smoothed_scores)*X),linewidth=1, color='g')
    plt.savefig(my_gene+'.png')
    print("Finished, written peak file and gene graph.")

if __name__ == "__main__":
    counts_bed, annot, N, X, min_gene_count, outfile_name, my_gene = main()
    if len(my_gene)>0:
        getSingleGenePeaks(counts_bed, annot, N, X, min_gene_count, outfile_name, my_gene)
    else:
        getAllPeaks(counts_bed, annot, N, X, min_gene_count, outfile_name)

