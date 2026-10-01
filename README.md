# panscape

panscape is a single Rust binary covering a pangenomics workflow — from raw
Nanopore/PacBio reads through alignment, filtering, motif search, VCF and
PAF comparison, GFA graph summarization, and genome annotation.

```
                                                            
  _ __     __ _   _ __    ___    ___    __ _   _ __     ___ 
 | '_ \   / _` | | '_ \  / __|  / __|  / _` | | '_ \   / _ \
 | |_) | | (_| | | | | | \__ \ | (__  | (_| | | |_) | |  __/
 | .__/   \__,_| |_| |_| |___/  \___|  \__,_| | .__/   \___|
 |_|                                          |_|           

panscape: analyzing pangenomes from reads to stats
       ************************************************
       Gaurav Sablok,
       gsablok@proton.me
      ************************************************

Usage: panscape <COMMAND>

Commands:
  fasta-convert             convert into fasta
  clipper-align             clipping the regions from the fastq
  scanner                   scans the reads for the motifs single occurence
  motifcatcher              motif plus upstream and the downstream
  selectedreads             selected reads writer
  filterreads               filter the reads prior to the length
  clip-seq                  remove the clip regions from the reads
  multi-clip-seq            remove the multitags for the fastqfile
  pangenome                 assemble pangenome
  minimap                   annotate reads
  stat                      annotated stats for your file
  pangenome-summarize       pangenome pre-computed alignment
  read-multisearch          multisearch reads across the reads
  harmonicmean              estimate the harmonic mean from the pangenome
  pangenome-matcher         pangenome matcher
  pan-arc                   pangenome annotator
  snatcher                  extract specific region from paf alignment
  precomputed-paf           generate stats from precomputed paf
  precompute-cds            extract the coding regions from the precomputed pangenome
  graph                     graph analyzer
  pangenome-bed             Pangenome bed constructor
  intergenic-noncoding      Intergenic extractor
  pan-reads-database        pangenome database
  bedtool-ancestral         analyze pangenome from the bedtools alignment to ancestral state
  vcf-aanalyze              analyze pangenome vcffiles
  pangenome-single-merge    merge single pangenome bedfile
  multi-bedtools-ancestral  bedtools sncestral multi-pangenome
  pafnet                    convert a paf alignment into a network/graph format (pafnet)
  help                      Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

## Build

```
cargo build --release
```

The binary is written to `target/release/panscape`. Every subcommand takes
a `thread` argument (a string, e.g. `"4"`) controlling how many threads its
internal rayon pool uses.

# panscape

panscape is a single Rust binary covering a pangenomics workflow — from raw
Nanopore/PacBio reads through alignment, filtering, motif search, VCF and
PAF comparison, GFA graph summarization, and genome annotation.

```
 _ __     __ _   _ __    ___    ___    __ _   _ __     ___
| '_ \   / _` | | '_ \  / __|  / __|  / _` | | '_ \   / _ \
| |_) | | (_| | | | | | \__ \ | (__  | (_| | | |_) | |  __/
| .__/   \__,_| |_| |_| |___/  \___|  \__,_| | .__/   \___|
|_|                                          |_|
```

## Build

```
cargo build --release
```

The binary is written to `target/release/panscape`. Every subcommand takes
a `thread` argument (a string, e.g. `"4"`) controlling how many threads its
internal rayon pool uses.

> If you're building with an older toolchain (rustc < 1.80), see
> [`OPTIMIZATIONS.md`](OPTIMIZATIONS.md) for a note on dependency pins.

## Commands and sample usage

Run `panscape --help` for the full list, or `panscape <command> --help` for
a single command's argument names. Sample commands below use placeholder
filenames — swap in your own paths.

### Reads / FASTQ processing

```bash
# Convert a FASTQ file to FASTA (writes fastaconvert.fasta)
panscape fasta-convert reads.fastq 4

# Report read-length distribution buckets (<=50kb, <=100kb, <=150kb, <=200kb, >200kb)
panscape stat reads.fastq 4

# Filter reads below a minimum length (writes filtered-ones.fasta)
panscape filterreads reads.fastq 5000 4

# Find every occurrence of a motif in each read (writes motipattern.txt)
panscape scanner reads.fastq ACGTACGT 4

# Find every regex match in each read (writes multi-searchpattern.txt)
panscape read-multisearch reads.fastq "ACG[TA]CGT" 4

# Grab a motif plus flanking upstream/downstream context
panscape motifcatcher reads.fastq ACGTACGT 50 50 4

# Keep only reads whose IDs appear in a list file (one ID per line)
panscape selectedreads reads.fastq wanted-ids.txt 4

# Remove a specific subsequence from each read
panscape clip-seq reads.fastq ACGTACGT 4

# Remove a region defined in a bed/region file from each read
panscape clipper-align reads.fastq region.bed 4

# Sequentially clip multiple tags listed in a text file, one per line
panscape multi-clip-seq reads.fastq multi-tags.txt 4
```

### Alignment / annotation

```bash
# Align reads against a protein set with minimap2 and slice out mRNA regions
panscape minimap reads.fastq proteins.fasta /usr/local/bin/minimap2 4

# Annotate a genome from a GFF + FASTA pair (writes mRNA.fasta, cds.fasta,
# cds-positive.fasta, cds-negative.fasta)
panscape pan-arc annotations.gff genome.fasta 4

# Extract CDS regions from a precomputed miniprot alignment
panscape precompute-cds alignment.paf reads.fasta 4

# Extract intergenic regions (see OPTIMIZATIONS.md — currently always empty)
panscape intergenic-noncoding alignment.gff reads.fasta 4
```

### PAF / pangenome comparison

```bash
# Summarize a single PAF alignment (aligned length, residue matches, etc.)
panscape pangenome-summarize alignment.paf 4

# Compare two PAF alignments and report matched query/reference pairs
panscape pangenome-matcher genome1.paf genome2.paf 4

# Slice matching regions out of query/reference FASTA using a PAF alignment
panscape snatcher alignment.paf query.fasta reference.fasta 4

# Summarize a precomputed PAF (query/reference totals, mean alignment length)
panscape precomputed-paf alignment.paf 4

# Harmonic mean estimate across per-line coordinate sets
panscape harmonicmean coordinates.txt 4

# Convert a PAF alignment into network/graph output (petgraph-backed,
# ported from pafnet-fast). mode is one of: net, edgelist, nodelist,
# rewrite, gexf. colors and prefixchar are only used by gexf mode (pass
# "" "" for the others).

# Pajek .net format (writes pafnet.net)
panscape pafnet alignment.paf net "" "" 4

# Plain edge list: src dst matches block_len mapq (writes pafnet.edges)
panscape pafnet alignment.paf edgelist "" "" 4

# id -> sequence name mapping (writes pafnet.nodes)
panscape pafnet alignment.paf nodelist "" "" 4

# Rewrite the PAF, substituting internal integer ids for names (writes pafnet.ids.paf)
panscape pafnet alignment.paf rewrite "" "" 4

# GEXF for Gephi, grouped/colored by name prefix (writes pafnet.gexf).
# colors.rgb is "group_name R G B" lines, matched against each sequence
# name's prefix split on the prefixchar (e.g. "group3.XXXXXX" with ".").
panscape pafnet alignment.paf gexf colors.rgb . 4
```

### GFA graph / pangenome BED

```bash
# Summarize a GFA graph (segment/link counts, total length, rank info)
panscape graph pangenome.gfa 4

# Build a BED-style summary from a GFA graph
panscape pangenome-bed pangenome.gfa 4

# Extract ancestral-state regions using bedtools + a reference FASTA
panscape bedtool-ancestral alignment.bed reference.fasta 100 /usr/local/bin/prank 4

# Same, for multiple alignments merged together
panscape multi-bedtools-ancestral alignment.bed reference.fasta /usr/local/bin/prank 4

# Merge a single pangenome BED file against a reference FASTA
panscape pangenome-single-merge alignment.bed reference.fasta 4
```

### VCF comparison

```bash
# Compare two VCF-derived range files against a reference FASTA
panscape vcf-aanalyze variants1.vcf variants2.vcf reference.fasta 4
```

### Assembly / database

```bash
# Run the hifiasm -> compleasm -> miniprot assembly pipeline (spawns external tools)
panscape pangenome hifi-reads.fastq 8 proteins.fasta

# Load reads into a local pangenome reads SQLite database
panscape pan-reads-database reads.fasta 4
```

## Verified against a real sample-files set

Every command below was actually run against a bundled `sample-files/`
directory (GFA graphs, PAF alignments, FASTA pairs, and a real PacBio HiFi
FASTQ read) and, where a reference output file shipped alongside the input,
diffed byte-for-byte against it.

```bash
# Graph summary — matches the shipped summarized_pangenome.txt exactly
$ panscape graph sample-files/pangenome-graph-analyze-sample-files/sample-pangenome.gfa 2
number_segment:8
number_links:11
number_arc:22
number_rank:2
total_segment_length:17572

# PAF comparison across two alignment files — output is byte-identical
# (matching md5sum) to the shipped comparative_query.txt /
# comparative_ref_write.txt reference files
$ panscape pangenome-matcher sample-files/pangenome-matcher-sample-files/test1.paf \
    sample-files/pangenome-matcher-sample-files/test2.paf 2

# Region extraction from a PAF + paired FASTA files — output is
# byte-identical to the shipped query-aligned.fasta / reference-aligned.fasta
$ panscape snatcher sample-files/pangenome-snatcher-sample-files/pangenome.paf \
    sample-files/pangenome-snatcher-sample-files/pangenome_query.fasta \
    sample-files/pangenome-snatcher-sample-files/pangenome_ref.fasta 2

# Real PacBio HiFi FASTQ (2 records, one with a genuine quality line, one
# truncated with no quality line at all) — both are parsed correctly with
# their own sequences, not mixed up:
$ panscape stat sample-files/samplepacbiohifi.fastq 4
Read length greater than 50KB:2	...

$ panscape scanner sample-files/samplepacbiohifi.fastq AAAAAAAAAA 4
$ cat motipattern.txt
"ERR10930361.1 magdelm64071_201030_115446/27"	34	44
"ERR10930361.1 magdelm64071_201030_115446/28"	34	44
```

Smaller synthetic checks for the FASTA-processing and motif-search commands:

```bash
$ cat test.fasta
@readA
AAAAACCCCCGGGGGTTTTT
@readB
CCCCCGGGGGTTTTTAAAAA
@readC
GGGGGTTTTTAAAAACCCCC

$ panscape scanner test.fasta CCCCC 2
$ cat motipattern.txt
"readA" 5   10
"readB" 0   5
"readC" 15  20
```


# Benchmarks: original vs. optimized, 15MB / 50MB / 100MB

## Setup

- Two binaries built from the same toolchain (rustc 1.75, same
  MSRV-compat dependency pins), so the comparison isolates the code changes
  in this repo, not a toolchain difference.
- **Original**: `git stash` back to the untouched source before any of the
  work in `OPTIMIZATIONS.md`, same `[profile.release]`-less Cargo.toml
  (Cargo's default release profile).
- **Optimized**: this repo's current state (buffered I/O, `fastareturn`
  fix, rayon parallelization, algorithmic fixes).
- Synthetic input: random 20,000bp reads, generated at three sizes in both
  plain 2-line FASTA form and 4-line FASTQ form (with real random Phred-like
  quality strings) — FASTQ is the case the original `fastareturn` handled
  worst, since every quality line was silently treated as extra sequence
  data.
- Each cell is the best-of-3 wall-clock time (`date +%s.%N` around the
  process, output redirected to `/dev/null`).
- **`--thread 1` throughout.** The sandbox this was run in has a single
  CPU core (`nproc` = 1), so rayon's parallelism has no real hardware to
  exploit here — every number below reflects only the non-threading
  changes (buffered writers, the `fastareturn` rewrite, eliminated
  redundant scans, O(n²)→O(n) fixes). On real multi-core hardware, running
  with `--thread N` for N > 1 should improve on these further for the
  commands that were parallelized; that couldn't be measured here.

## Results

| format | size  | command       | original (s) | optimized (s) | speedup |
|--------|-------|---------------|---------------|----------------|---------|
| fasta  | 15MB  | filterreads   | 0.039         | 0.028          | 1.37x   |
| fasta  | 15MB  | scanner       | 0.033         | 0.014          | 2.32x   |
| fasta  | 15MB  | stat          | 0.023         | 0.014          | 1.67x   |
| fasta  | 50MB  | filterreads   | 0.117         | 0.101          | 1.16x   |
| fasta  | 50MB  | scanner       | 0.095         | 0.046          | 2.05x   |
| fasta  | 50MB  | stat          | 0.087         | 0.047          | 1.87x   |
| fasta  | 100MB | filterreads   | 0.246         | 0.207          | 1.18x   |
| fasta  | 100MB | scanner       | 0.175         | 0.093          | 1.87x   |
| fasta  | 100MB | stat          | 0.175         | 0.097          | 1.81x   |
| fastq  | 15MB  | filterreads   | 0.033         | 0.018          | 1.86x   |
| fastq  | 15MB  | scanner       | 0.029         | 0.014          | 2.09x   |
| fastq  | 15MB  | stat          | 0.029         | 0.013          | 2.14x   |
| fastq  | 50MB  | filterreads   | 0.064         | 0.060          | 1.05x   |
| fastq  | 50MB  | scanner       | 0.063         | 0.044          | 1.43x   |
| fastq  | 50MB  | stat          | 0.063         | 0.045          | 1.40x   |
| fastq  | 100MB | filterreads   | 0.121         | 0.128          | **0.94x** |
| fastq  | 100MB | scanner       | 0.121         | 0.095          | 1.26x   |
| fastq  | 100MB | stat          | 0.121         | 0.094          | 1.28x   |

## The one regression: `filterreads` on large FASTQ input

At 100MB FASTQ, `filterreads` is about 6% *slower* on the optimized build
(0.121s → 0.128s), and only a marginal 1.05x at 50MB. Every other
command/size/format combination improved. This was checked across 5 repeat
runs each (not a one-off noisy sample) and is a real, reproducible effect
in this single-core environment, not measurement noise.
