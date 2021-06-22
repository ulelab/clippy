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
```
user$ python clip.py -h
usage: clip.py [-h] -i INPUTBED -o OUTPUTPREFIX -a ANNOT [-n [WINDOWSIZE]]
               [-x [ADJUST]] [-mg [MINGENECOUNTS]] [-mb [MINPEAKCOUNTS]]
               [-d [DISTANCE]] [-g [MYGENE]]

Call CLIP peaks.

required arguments:
  -i INPUTBED, --inputbed INPUTBED
                        bed file containing crosslink counts at each position
  -o OUTPUTPREFIX, --outputprefix OUTPUTPREFIX
                        prefix for output files
  -a ANNOT, --annot ANNOT
                        gtf annotation file

optional arguments:
  -h, --help            show this help message and exit
  -n [WINDOWSIZE], --windowsize [WINDOWSIZE]
                        rolling mean window size [DEFAULT 50]
  -x [ADJUST], --adjust [ADJUST]
                        adjustment for prominence [DEFAULT 1]
  -hc [HEIGHT_CUTOFF], --height_cutoff [HEIGHT_CUTOFF]
                        proportion of prominence [DEFAULT 0.8]
  -mg [MINGENECOUNTS], --mingenecounts [MINGENECOUNTS]
                        min counts per gene to look for peaks [DEFAULT 5]
  -mb [MINPEAKCOUNTS], --minpeakcounts [MINPEAKCOUNTS]
                        min counts per broad peak [DEFAULT 5]
  -g [MYGENE], --mygene [MYGENE]
                        gene name, limits analysis to single gene
  -t [THREADS], --threads [THREADS]
                        number of threads to use
  -cf [CHUNKSIZE_FACTOR], --chunksize_factor [CHUNKSIZE_FACTOR]
                        A factor used to control the number of jobs given to a thread at a time. A larger number reduces the number of jobs per chunk. Only increase if you experience crashes [DEFAULT 4]
  -int, --interactive   starts a Dash server to allow for interactive parameter tuning
```
*A note on -g, my-gene option*

When you are testing parameters you might want to check what they look like on certain genes, to save a graph of a given gene provide the name or reference (as in your annotation) and the code will only run for your given gene and output a graph like the one below in the 'concept' section. It will save with the gene name and parameters as a file name.

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


### Concept
Using the annotation provided, crosslinks over each gene are smoothed using a rolling mean. The window can be decided by the user. For each gene the mean of the smoothed signal is taken (red line) and the mean + (standard deviation * adjustment factor) (green line) is taken. The mean is used to define the minimum height of a peak. The mean + (standard deviation * adjustment factor) is taken to define the minimum prominence of a peak. Please see [here](https://en.wikipedia.org/wiki/Topographic_prominence#:~:text=The%20prominence%20of%20a%20peak,or%20key%20saddle%2C%20or%20linking) for the definition of topographical prominence. Essentially this parameter allows that we do not call many shallow peaks in a region where there is a clearly more prominent peak. 

In the image below, you can see the four positions the algorithm picks out in this gene as peaks. To get broader windows we take these single nt positions and merge adjacent crosslinks.

![Image of gene](pmt2_demo.png)


### Developer Functions

If you plan to contribute to the Clippy code we have some helpful functions for development. To run the automated testing, use:

```
pytest --cov=clip
```

You might be interested in how long certain functions take to run. To run the profiling code:

```
pytest -k profiling --profile-svg
python -m gprof2dot -f pstats prof/get_the_peaks.out | dot -Tpdf -o prof/get_the_peaks.pdf
```


### Authors
Charlotte Capitanchik - charlotte.capitanchik@crick.ac.uk
Marc Jones - marc.jones@crick.ac.uk
