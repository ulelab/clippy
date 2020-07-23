# clippy
A wrapper around scipy "find_peaks" function to enable peak calling of CLIP data.

### Requirements
This is a python 3 program requiring the following modules:
 - numpy
 - scipy
 - pybedtools
 - pandas
 - matplotlib


### Usage
```
get_peaks.py -i <input-bed-counts> -o <output-file-prefix> -a <annotation-gtf> \
[optional... -n <rolling-mean-window> -x <adjust-prominence> -m <min-gene-count>]
```

### Concept
