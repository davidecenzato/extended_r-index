# Extended r-index

The extended r-index is a tool implementing the extended r-index data structure for string collections.
The code id adapted from the r-index by Nicola Prezza, Travis Gagie and Gonzalo Navarro (https://github.com/nicolaprezza/r-index). 

# Usage

### Construction of the extended r-index:
```
usage: ext_r-index.py [-h] [--construct] [-w WSIZE] [-p MOD] [-b B] [--first] [--pfile PFILE] [--count] [--locate] [--verbose] input

Tool to build the extended r-index of string collections.

positional arguments:
  input                 input fasta file name

options:
  -h, --help            show this help message and exit
  --construct           constructs the extended r-index (def. False)
  -w WSIZE, --wsize WSIZE
                        sliding window size for the PFP algorithm (def. 10)
  -p MOD, --mod MOD     hash modulus for PFP (def. 100)
  -b B, --B B           bitvector block size for predecessor queries (def. 2)
  --first               sample first rotation of each sequence (def. False)
  --pfile PFILE         pattern file path (def. <input filename.pat>)
  --count               compute count queries (def. False)
  --locate              compute locate queries (def. False)
  --verbose             verbose (def. False)
```
The extended r-index construction using the cyclic PFP algorithm is enabled using the `--construction` flag. The count and locate queries computation
is enabled using the `--count` and `--locate` flag, the file containing the patterns, in fasta format, is defined using the `--pfile` flag.

# Example
### Download and Compile

```console
git clone https://github.com/davidecenzato/extended_r-index.git
cd extended_r-index
mkdir build
cd build
cmake ..
make
```

### Run on Example Data

```console
// Build the extended r-index of a genomic string collection
python3 ext_r-index.py data/yeast.fasta --construct -w 10 -p 100  
// Generate a pattern file (100 patterns of length 100)
build/genpattern data/yeast.fasta 100 100 data/yeast.patt 1 1
// Run count queries
python3 ext_r-index.py data/yeast.fasta --count --pfile data/yeast.patt  
```
# External resources

* [pfpebwt](https://github.com/davidecenzato/PFP-eBWT.git)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
* [gSACA-K](https://github.com/felipelouza/gsa-is.git)
* [malloc_count](https://github.com/bingmann/malloc_count)

# Authors

### Theoretical results:

* Christina Boucher
* Davide Cenzato
* Zsuzsanna Lipták
* Massimiliano Rossi
* Marinella Sciortino

### Implementation and experiments:

* [Davide Cenzato](https://github.com/davidecenzato) 
* [Massimiliano Rossi](https://github.com/maxrossi91)