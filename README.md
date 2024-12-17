Rust program to normalize coverage between pairs of sequence regions (e.g. genes) in BAM files

This requires at least two bam files (sorted & indexed) and a tab separated list of regions (see example).  There must be one region per line per bam file.

Normalizes coverage between regions given by a tab separated file
ie chr:1-1000<tab>chr2:1000-2000<tab>...

The minimum coverage setting will cause regions which fall below the threshold to be skipped and the matching regions will instead be normalized to the lowest which falls above the threshold. For example:

--min_covg 5 on the following<br>
Reg1: 7x, Reg2: 3x, Reg3: 10x, Reg4: 8x<br>
would yield<br>
Reg1: 7x, Reg2: 3x, Reg3: 7x, Reg4: 7x

Usage: normalize_paired_regions [OPTIONS] --regions <PATH> --bams <PATHS>...

Options:
  -r, --regions <PATH>       Path to tab-separated regions file
  -b, --bams <PATHS>...      Paths to bam files, must match # of cols in regions file! ./bam_1 ./bam_2 [... ./bam_n]
  -u, --allow_unpaired       Allows single reads through
  -v, --verbose              prints debugging info
  -c, --min_covg <min_covg>  Minimum coverage (default 0)
  -h, --help                 Print help
