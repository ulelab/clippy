import numpy as np
import scipy.signal as sig
from scipy.ndimage.filters import uniform_filter1d
from scipy import stats
import pybedtools
import pandas as pd
import matplotlib.pyplot as plt
import getopt, sys

def main(argv):
   counts_bed = ''
   annot = ''
   N = 50 # rolling mean smooth
   X = 1 # multiplication factor for standard deviation - to define prominence
   min_gene_count = 5 # minimum counts on a gene to consider it for peak calling
   outfile_prefix = ''
   my_gene = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:a:nxmg",["counts_bed=","outfile_prefix=","annot=","N=","X=","min_gene_count=","my_gene"])
   except getopt.GetoptError:
      print 'get_peaks.py -i <input-bed-counts> -o <output-file-prefix> -a <annotation-gtf> [optional... -n <rolling-mean-window> -x <adjust-prominence> -m <min-gene-count>]'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'get_peaks.py -i <input-bed-counts> -o <output-file-prefix> -a <annotation-gtf> [optional... -n <rolling-mean-window> -x <adjust-prominence> -m <min-gene-count>]'
         sys.exit()
      elif opt in ("-i", "--counts_bed"):
         counts_bed = arg
      elif opt in ("-o", "--outfile_prefix"):
         outfile_prefix = arg
      elif opt in ("-a", "--annot"):
         annot = arg
      elif opt in ("-n", "--N"):
         N = arg
      elif opt in ("-x", "--X"):
         X = arg
      elif opt in ("-m", "--min_gene_count"):
         min_gene_count = arg
      elif opt in ("-g", "--my_gene"):
         my_gene = arg
      outfile_name=outfile_prefix+str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".bed"
   return(counts_bed, annot, int(N), int(X), int(min_gene_count), outfile_name, my_gene)

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
    return(peaks_in_gene, roll_mean_smoothed_scores)

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
    ang = ang[ang.attributes.str.upper.contains(my_gene.upper())]

    if len(ang) == 0:
        sys.exit("Couldn't find your gene in the provided annotation")
    elif len(ang) > 1:
        sys.exit("Found more than one gene containing that name - could you be more specific?")

    ang = pybedtools.BedTool.from_dataframe(annot_gene)
 
    goverlaps = pho92_iclip.intersect(ang, s=True, wo=True).to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand','chrom2','source','feature','gene_start', 'gene_stop','nothing','strand2','nothing2','gene_name','interval'])
    goverlaps.drop(['name','chrom2','nothing','nothing2','interval','strand2','source','feature'], axis=1, inplace=True)

    peaks, roll_mean_smoothed_scores = getThePeaks(goverlaps)

    if peaks.empty:
        sys.exit("No peaks found in this gene with the current parameters")

    outfile_name=my_gene+str(N)+"_stdev"+str(X)+"_minGeneCount"+str(min_gene_count)+".bed"
    peaks.to_csv(outfile_name,sep="\t",header=False,index=False)

    # Make graph of gene
    plt.plot(roll_mean_smoothed_scores, '-bD', markevery=peaks[0].tolist())
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