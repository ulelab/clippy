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


### Usage
```
clip.py -i <input-bed-counts> -o <output-file-prefix> -a <annotation-gff> \
[optional... -n <rolling-mean-window> (default: 50) \
-x <adjust-prominence> (default: 1) \
-m <min-gene-count> (default: 5) \
-g <my-gene> (default: NULL)]
```
*A note on -g, my-gene option*

When you are testing parameters you might want to check what they look like on certain genes, to save a graph of a given gene provide the name or reference (as in your annotation) and the code will only run for your given gene and output a graph like the one below in the 'concept' section. It will save with the gene name and parameters as a file name.

*A note on annotation gff*

The code only requires that you have a feature labelled "gene" in the 3rd column of your gff, and assumes that the 9th column of your gff will uniquely identify your genes and contain some kind of gene name or ID.

### Run test data

```
clip.py -i test/crosslinkcounts.bed -o test -a test/annot.gff \
-n 50 -x 1 -m 5 -g pmt2
```

### Concept
Using the annotation provided, crosslinks over each gene are smoothed using a rolling mean. The window can be decided by the user. For each gene the mean of the smoothed signal is taken (red line) and the mean + (standard deviation * adjustment factor) (green line) is taken. The mean is used to define the minimum height of a peak. The mean + (standard deviation * adjustment factor) is taken to define the minimum prominence of a peak. Please see [here](https://en.wikipedia.org/wiki/Topographic_prominence#:~:text=The%20prominence%20of%20a%20peak,or%20key%20saddle%2C%20or%20linking) for the definition of topographical prominence. Essentially this parameter allows that we do not call many shallow peaks in a region where there is a clearly more prominent peak. 

In the image below, you can see the four positions the algorithm picks out in this gene as peaks. To get broader windows we take these single nt positions and merge adjacent crosslinks.

![Image of gene](pmt2_demo.png)


### Development ideas
- [ ] Check for overlapping genes and have a heirarchy to deal with them.
- [ ] Try instead of merging crosslinks to get broader peak regions, use the intercept of the smoothed crosslink peaks with the green prominence line. Note: I have tried using the intrinsic width property of "find_peaks" and this seems to call things way too wide.
- [ ] Parallelise looping through genes to improve speed.

### Author
Charlotte Capitanchik - charlotte.capitanchik@crick.ac.uk
