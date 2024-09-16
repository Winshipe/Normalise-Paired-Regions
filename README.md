Rust program to normalize coverage between pairs of sequence regions (e.g. genes) in BAM files

Uses rust_htslib (https://docs.rs/rust-htslib/latest/rust_htslib/index.html) and CLAP

This requires two bam files (sorted & indexed) and a tab separated list of regions (see example)

Usage: normalize_paired_regions [OPTIONS] --bam_1 --bam_2 --regions

Options:

-1, --bam_1 Path to first BAM file

-2, --bam_2 Path to second BAM file

-r, --regions Path to tab-separated regions file

-o, --output_1 Output path for 1st BAM file (optional)

-n, --output_2 Output path for 1st BAM file (optional)

-c, --min_covg <min_covg> Minimum coverage (default 0)

-h, --help Print help
