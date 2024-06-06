# Sperm Seq Scripts

## Name
Sperm-seq scripts

## Description
Some helper scripts for manipulating Sperm-seq data from [Bell, et al. 2020 "Insights into variation in meiosis from 31,228 human sperm genomes"](https://doi.org/10.1038/s41586-020-2347-0)

**download_and_zip.sh**
- sequentially submits jobs to SLURM to download, then zip *.fastq files from prefetched *.sra files
- requires list of *.sra files along with .ngc file. More info [here](https://www.ncbi.nlm.nih.gov/sra/docs/sra-dbgap-download/)

**enumerate_nuisance_barcodes.pl**
- extracts common contaminants from fastqc's [contaminant_list.txt](https://github.com/s-andrews/FastQC/blob/master/Configuration/contaminant_list.txt)
- enumerates all 14mers <= 1 edit from common contaminants, their reverses, complements, and reverse complements
- enumerates all 14mers <= 1 edit from homopolymer sequences (e.g `GGGGGGGGGGGGGG` and `GGGGTGGGGGGGGG`)
- use output as a blocklist for `fix_bam_tags.pl`

**fix_bam_tags.pl**
- extracts the cell barcode from the XM tag and the umi from the XC tag
- creates new RG tags based on the donor id, run id, and cell barcode
- creates CB tag with cell barcode sequence
- creates MI and OX tags with the umi sequence
- filters low abundance cell barcodes (default <= 1)
- filters barcodes according to blocklist
- filters out reads marked as duplicates
- produces new sorted bam file
- counts reads aligning to CHR1/PAR1/NPX/PAR2/NPY
- calculates RPKM for each region/readgroup
- calls each sample as:
  * X-bearing (1) or Y-bearing (0)
  * Sex chromosome aneuploid (1) or euploid (0)
  * Doublet (1) or singlet (0)
- output TSV with calls and RPKMs for each RG for use as karyotype calls input to `identify_crossovers.pl`

**filter_x_bearing.pl**

- processes TSV of karyotype calls from `fix_bam_tags.pl` on stdin
- prints only those lines with
  * chrX coverage >0.5x and < 2x chr1 coverage
  * chrY coverage <0.5x chrX coverage
  * chrY coverage <0.25x chr1 coverage

**filter_y_bearing.pl**

- processes TSV of karyotype calls from `fix_bam_tags.pl` on stdin
- prints only those lines with
  * chrY coverage >0.5x and < 2x chr1 coverage
  * chrX coverage <0.5x chrX coverage
  * chrX coverage <0.25x chr1 coverage

**pipeline.pl**

- processes SRA files from SraRunTable.csv one at a time (by donor):
  * convert SRA to paired fastq.gz
  * convert paired fastq.gz to unaligned bam
  * add cell barcode and UMI tags to unaligned bam
  * extract single-ended fastq | cutadapt | bowtie2 > sorted bam
  * merge unaligned and aligned bams
  * mark duplicate reads (using UMIs)
  * fix bam tags
  * for each run: save fixed bam, tsv of karyotype and aneuploidy calls, mkdup stats
  * for each donor: merge fixed bams for all runs, run freebayes to generate vcf

**identify_crossovers.pl**
- reads in:
  * PAR VCF
  * table of sex-chromosome karyotype calls per sample (X-bearing or not)
  * a list of samples to plot
  * chromosome name to phase and plot (e.g. NC_000023.11)
- filters variants for:
  * quality > 30
  * SNV only
  * At least two REF and ALT calls
  * 0.25 < observed allele frequency < 0.75
  * (optionally) a list of banned unreliable markers
  * (optionally) consistency in NCO samples
- reconstructs donor haplotypes in 2 passes:
  * Using list of known (or inferred) X- and Y-bearing sperm
  * starting at PAB, iteratively phase sites toward TEL
  * assign alleles to haplotypes by minimizing number of sperm that must change haplotype from previous sites
- filters out:
  * Double crossovers supported by a single SNP (likely genotyping error)
  * Tightly spaced (~5 cM) double crossovers (possible gene conversion)
- prints a bedGraph of recombination rate (cM/Mb)
- prints a SVG image of schematic chromosomes with colored tickmarks for phased variants
- prints a TSV file with crossover locations

**average_bedgraph.pl**

- supply a list of bedGraphs files as arguments
- prints a new bedGraph, averaging the values for every base

**plot_par_recombination.py**
- requirements:
  * python 3.10.9
  * numpy 1.26.4
  * pandas 2.2.1
  * matplotlib 3.8.0
- reads in:
  * parent directory to find all recombination rates .txt files
- generates:
  * plot of average recombination rate across the human PAR1
  * plot of recombination rate across human PAR1 for each donor
  * histogram plot of recombination AUC for each donor for:
    - entire PAR1
    - interval between proposed and established PAR1 boundaries
- outputs:
  * calculated mean AUC and standard deviation across:
    - entire PAR1
    - interval between proposed and established PAR1 boundaries

## Installation
Download these perlscripts and put them in your path

## Usage

**fix_bam_tags.pl**
```
Usage: ./fix_bam_tags.pl -i input.bam -o output.bam -t rpkm.tsv

Mandatory Arguments:

  -i  bam file
  -o  sorted bam file
  -t  tsv of reads per kilobase per million mapped

Optional Arguments:

  -u  print usage information and exit
  -h  print this help and exit

  For customizing read group id:
  -d  donor id
  -r  run id

  Filtering cell barcodes:
  -m  minimum counts
  -b  blocklist (one barcode per line)
```

**pipeline.pl**

```
Usage: ./pipeline.pl -c SraRunTable.csv -s /directory/with/sra_files -b /path/to/bowtie2/index -g genome.fasta -n nuisance.txt -v dbsnp.vcf -o /directory/for/output

Mandatory Arguments:

  -c  SraRunTable.csv from ncbi
  -s  directory with sra files
  -b  path to bowtie2 index
  -g  fasta file of genome
  -n  text file containing nuisance cell barcodes (1 per line)
  -v  vcf of SNPs to call
  -o  directory for output

Optional Arguments:

  -f  directory with fastq.gz files
  -u  print usage information and exit
  -h  print this help and exit
```

**identify_crossovers.pl**
```
Usage: identify_crossovers.pl -v input.vcf -t calls.tsv -s samples.txt -c chrX -p 2781480

Mandatory Arguments:

  -v  vcf file
  -k  tsv with sex chromosome karyotype calls
  -s  list of samples (one per line)
  -c  sex chromosome name (must match VCF)
  -p  NPX/NPY position (integer)

Optional Arguments:

  -u  print usage information and exit
  -h  print this help and exit

  Outputs:

  -i name of image file (.svg)
  -l table of crossover locations (.tsv)
  -g bedGraph output file (.txt)
```

**plot_par_recombination.py**
```
usage: ./plot_averaged_par_recomb.py [-h] [-d DIR] [-o OUT] [-c COLS] [-w WIDTH] [-l LENGTH]

Find recombination files and plot results in png and svg form.

optional arguments:
  -h, --help            show this help message and exit
  -d DIR, --dir DIR     Path to parent directory to grab files. Defaults to current folder.
  -o OUT, --out OUT     Path to directory to save. Defaults to current folder.
  -c COLS, --cols COLS  Number of columns for donor separated graphs. Defaults to 2.
  -w WIDTH, --width WIDTH
                        Figure width of donor separated graphs. Defaults (8) inches.
  -l LENGTH, --length LENGTH
                        Figure length (height) of donor separated graphs. Defaults (10) inches.
  -f FONT, --font FONT  Path to .ttf font file (set up for Helvetica). Defaults to DejaVu Sans.
```

## Support
Leave an issue here.

## License
MIT License
