# clippy
A wrapper around scipy "[find_peaks](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html)" function to enable peak calling of CLIP data.
![A dumb joke](smallerclippy.png)

### Requirements
This is a python 3 program requiring the following modules:
 - numpy
 - scipy
 - pybedtools
 - pandas
 - matplotlib

It is important to use bedtools v2.26.0. Also note, running matplotlib on a cluster it is helpful to add `export QT_QPA_PLATFORM='offscreen'` to your bash profile to avoid display errors.

If you use conda, we provide an environment.yaml:
 
```
conda env create -f environment.yml
conda activate clippy
```

### Usage

Clippy is now available to install from Bioconda:
```
conda install -c bioconda clippy
```

```
usage: clip.py [-h] [-v] -i INPUTBED -o OUTPUTPREFIX -a ANNOT -g GENOME_FILE
               [-n [WINDOWSIZE]] [-up [UPSTREAM_EXTENSION]]
               [-down [DOWNSTREAM_EXTENSION]] [-x [ADJUST]]
               [-hc [HEIGHT_CUTOFF]] [-mg [MINGENECOUNTS]]
               [-mb [MINPEAKCOUNTS]] [-m [MYGENE]] [-t [THREADS]]
               [-cf [CHUNKSIZE_FACTOR]] [-int] [-nei]
               [-inter [INTERGENIC_PEAK_THRESHOLD]] [-alt [ALT_FEATURES]]

Call CLIP peaks.

required arguments:
  -i INPUTBED, --inputbed INPUTBED
                        bed file containing crosslink counts at each position
  -o OUTPUTPREFIX, --outputprefix OUTPUTPREFIX
                        prefix for output files
  -a ANNOT, --annot ANNOT
                        gtf annotation file
  -g GENOME_FILE, --genome_file GENOME_FILE
                        genome file containing chromosome lengths, used by
                        BEDTools for genomic operations

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -n [WINDOWSIZE], --windowsize [WINDOWSIZE]
                        rolling mean window size [DEFAULT 50]
  -up [UPSTREAM_EXTENSION], --upstream_extension [UPSTREAM_EXTENSION]
                        upstream extension added to gene models [DEFAULT 0]
  -down [DOWNSTREAM_EXTENSION], --downstream_extension [DOWNSTREAM_EXTENSION]
                        downstream extension added to gene models [DEFAULT 0]
  -x [ADJUST], --adjust [ADJUST]
                        adjustment for prominence [DEFAULT 1]
  -hc [HEIGHT_CUTOFF], --height_cutoff [HEIGHT_CUTOFF]
                        proportion of prominence [DEFAULT 0.8]
  -mg [MINGENECOUNTS], --mingenecounts [MINGENECOUNTS]
                        min counts per gene to look for peaks [DEFAULT 5]
  -mb [MINPEAKCOUNTS], --minpeakcounts [MINPEAKCOUNTS]
                        min counts per broad peak [DEFAULT 5]
  -m [MYGENE], --mygene [MYGENE]
                        gene name, limits analysis to single gene
  -t [THREADS], --threads [THREADS]
                        number of threads to use
  -cf [CHUNKSIZE_FACTOR], --chunksize_factor [CHUNKSIZE_FACTOR]
                        A factor used to control the number of jobs given to a
                        thread at a time. A larger number reduces the number
                        of jobs per chunk. Only increase if you experience
                        crashes [DEFAULT 4]
  -int, --interactive   starts a Dash server to allow for interactive
                        parameter tuning
  -nei, --no_exon_info  Turn off individual exon and intron thresholds
  -inter [INTERGENIC_PEAK_THRESHOLD], --intergenic_peak_threshold [INTERGENIC_PEAK_THRESHOLD]
                        Intergenic peaks are called by first creating
                        intergenic regions and calling peaks on the regions as
                        though they were genes. The regions are made by
                        expanding intergenic crosslinks and merging the
                        result. This parameter is the threshold number of
                        crosslinks required to include a region. If set to
                        zero (default), no intergenic peaks will be called.
                        When using this mode, the intergenic regions used will
                        be output as a GTF file. [DEFAULT 0]
  -alt [ALT_FEATURES], --alt_features [ALT_FEATURES]
                        A list of alternative GTF features to set individual
                        height thresholds on in the comma-separated format
                        <alt_feature_name>-<gtf_key>-<search_pattern>
```

*A note on annotation gff*

The code only requires that you have a feature labelled "gene" in the 3rd column of your gff, and assumes that the 9th column of your gff will uniquely identify your genes and contain some kind of gene name or ID.

### Run test data

Get peaks for one gene along with an image of that gene.

```
python clip.py -i tests/data/crosslinkcounts.bed -o test_plot -a tests/data/annot.gff \
-n 50 -x 1 -mg 5 -mb 5 -g pmt2 -hc 0.8
```

Start a test instance of the interactive parameter search server:

```
python clip.py -i tests/data/crosslinkcounts.bed -o TESTING -a tests/data/annot.gff -int
```

#### Test data generation

For CEP295 data:

```
cat tests/data/gencode.v38.annotation.gtf | awk '{if($1=="chr11" && 93661682<=$4 && $5<=93730358){print($0)}}' > tests/data/gencode.v38.cep295.gtf

gunzip -c tests/data/tardbp-egfp-hek293-2-20201021-ju_mapped_to_genome_reads_single.bed.gz | awk '{if($1=="chr11" && 93661682<=$2 && $3<=93730358){print($0)}}' > tests/data/cep295.bed
```

For RBFOX2 data:

```
cat tests/data/gencode.v38.annotation.gtf | awk '{if($1=="chr19" && 10000000<=$4 && $5<=11000000){print($0)}}' > tests/data/gencode.v38.chr19_10M_11M.gtf

gunzip -c tests/data/HepG2_RBFOX2.xl.bed.gz | awk '{if($1=="chr19" && 10000000<=$2 && $3<=11000000){print($0)}}' > tests/data/rbfox2_chr19_10M_11M.bed
```

### Concept
Using the annotation provided, crosslinks over each gene are smoothed using a rolling mean. The window can be decided by the user. For each gene the mean of the smoothed signal is taken (red line) and the mean + (standard deviation * adjustment factor) (green line) is taken. The mean is used to define the minimum height of a peak. The mean + (standard deviation * adjustment factor) is taken to define the minimum prominence of a peak. Please see [here](https://en.wikipedia.org/wiki/Topographic_prominence#:~:text=The%20prominence%20of%20a%20peak,or%20key%20saddle%2C%20or%20linking) for the definition of topographical prominence. Essentially this parameter allows that we do not call many shallow peaks in a region where there is a clearly more prominent peak. 

In the image below, you can see the four positions the algorithm picks out in this gene as peaks. To get broader windows we take these single nt positions and merge adjacent crosslinks.

![Image of gene](pmt2_demo.png)


### Developer Functions

If you plan to contribute to the Clippy code we have some helpful functions for development. To run the automated testing, use:

```
pytest --cov=clip -k "not profiling and not web"
```

You might be interested in how long certain functions take to run. To run the profiling code:

```
pytest -k profiling --profile-svg
python -m gprof2dot -f pstats prof/get_the_peaks.out | dot -Tpdf -o prof/get_the_peaks.pdf
```


### Authors
Charlotte Capitanchik - charlotte.capitanchik@crick.ac.uk
Marc Jones - marc.jones@crick.ac.uk
