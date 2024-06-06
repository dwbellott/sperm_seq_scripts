#!/usr/bin/env perl

use strict;
use Cwd;
use Cwd 'abs_path';
use File::Which;
use File::Path 'rmtree';
use File::Find;
use Text::CSV qw( csv );
use Getopt::Long;
use vars qw/$VERSION/;

BEGIN {
  $VERSION = '0.1'
};

sub main ();
sub usage ($);
sub help ();
sub find_executables ();
sub process_run ($$$);

&main ();

sub main (){


  #setup command-line options
  my (
    $csv,
    $sra,
    $bowtie,
    $genome,
    $nuisance,
    $vcf,
    $output,
    $fastq,
    $help,
    $usage);

  if (!GetOptions(
      'csv=s' => \$csv,
      'sra=s' => \$sra,
      'bowtie=s' => \$bowtie,
      'genome=s' => \$genome,
      'nuisance=s' => \$nuisance,
      'vcf=s' => \$vcf,
      'output=s' => \$output,
      'fastq=s' => \$fastq,
      'help'  => \$help,
      'usage' => \$usage)){
    print STDERR &usage("$0 (Version: $VERSION)");
    exit -1;
  }

  if ($help){
    print STDOUT &help;
    exit;
  }

  if ($usage){
    print STDOUT &usage("(Version: $VERSION)");
    exit;
  }

  if (defined($csv)){
    if (-e $csv && -f _ && -r _){
      $csv = abs_path($csv);
    }else{
      print STDOUT &usage("SraRunTable.csv at $csv is not readable\n");
      exit;
    }
  }else{
    print STDOUT &usage("Please supply the location of the SraRunTable.csv with -c\n");
    exit;
  }

  if (defined($sra)){
    if (-e $sra && -d _ && -r _){
      $sra = abs_path($sra);
    }else{
      print STDOUT &usage("directory at $sra is not readable\n");
      exit;
    }
  }else{
    print STDOUT &usage("Please supply a directory containing sra files with -s\n");
    exit;
  }

  if (defined($bowtie)){
    if (-e $bowtie.".1.bt2" && -f _ && -r _ &&
        -e $bowtie.".2.bt2" && -f _ && -r _ &&
        -e $bowtie.".3.bt2" && -f _ && -r _ &&
        -e $bowtie.".4.bt2" && -f _ && -r _ &&
        -e $bowtie.".rev.1.bt2" && -f _ && -r _ &&
        -e $bowtie.".rev.2.bt2" && -f _ && -r _ ){
      $bowtie = abs_path($bowtie.".1.bt2");
      $bowtie =~ s/\.1\.bt2$//;
    }else{
      print STDOUT &usage("bowtie2 index at $bowtie is not readable\n");
      exit;
    }
  }else{
    print STDOUT &usage("Please supply a bowtie2 index with -b\n");
    exit;
  }

  if (defined($genome)){
    if (-e $genome && -f _ && -r _ ){
      $genome = abs_path($genome);
    }else{
      print STDOUT &usage("genome fasta at $genome is not readable\n");
      exit;
    }
  }else{
    print STDOUT &usage("Please supply a genome fasta with -g\n");
    exit;
  }

  if (defined($vcf)){
    if (-e $vcf && -f _ && -r _ ){
      $vcf = abs_path($vcf);
    }else{
      print STDOUT &usage("vcf file at $vcf is not readable\n");
      exit;
    }
  }else{
    print STDOUT &usage("Please supply a vcf file with trusted SNPs using -v\n");
    exit;
  }

  if (defined($output)){
    if (-e $output){
      if (-d $output && -w _){
        $output = abs_path($output);
      }else{
        print STDOUT &usage("directory at $output is not writable\n");
        exit;
      }
    }else{
      $output = abs_path($output);
      mkdir($output);
    }
  }else{
    print STDOUT &usage("Please supply an output directory with -o\n");
    exit;
  }

  if (defined($nuisance)){
    if (-e $nuisance && -f _ && -r _){
      $nuisance = abs_path($nuisance);
    }else{
      print STDOUT &usage("File of nuisance barcodes at $nuisance is not readable\n");
      exit;
    }
  }else{
    print STDOUT &usage("Please supply nuisance barcodes with -n\n");
    exit;
  }

  if (defined($fastq)){
    if (-e $fastq && -d _ && -r _){
      $fastq = abs_path($fastq);
    }else{
      print STDOUT &usage("directory at $fastq is not readable\n");
      exit;
    }
  }

  my $environment = &find_executables();

  $environment->{'SraRunTable'} = $csv;
  $environment->{'sra_directory'} = $sra;
  $environment->{'bowtie2_index'} = $bowtie;
  $environment->{'genome_fasta'} = $genome;
  $environment->{'db_snp_vcf'} = $vcf;
  $environment->{'nuisance_barcodes'} = $nuisance;
  $environment->{'output_directory'} = $output;
  $environment->{'fastq_directory'} = $fastq;

  foreach my $p (keys(%{$environment})){
    print join("\t", $p, $environment->{$p})."\n";
  }

  my $dict = $genome;
  $dict =~ s/\.[a-z]$/dict/;
  unless (-e $dict){
    my $cmd = $environment->{'picard'}." CreateSequenceDictionary R=".$genome." O=".$dict;
    print STDOUT "$cmd\n";
    system($cmd);
  }
  $environment->{'genome_dictionary'} = $dict;

  # set up some variables;
  my %donor;
  my %run;
  my %donor_by_run;
  my %runs_by_donor;
  # read in SraRunTable from NCBI:
  my $aoh = csv (in => $csv, headers => "auto");
  foreach my $h (@{$aoh}){
    my $r = $h->{'Run'};
    my $d = $h->{'submitted_subject_id'};
    $donor{$d} = 1;
    $run{$r} = 1;
    $donor_by_run{$r} = $d;
    push(@{$runs_by_donor{$d}}, $r);
  }

  # read in a table of information
  my @donors = keys(%donor);
  my @runs = sort {$a cmp $b} keys(%run);


  #processing steps for each donor:
  foreach my $d (@donors){
    my @bams;
    my $cmd;
    ## process each run to get a fixed BAM
    foreach my $r (@{$runs_by_donor{$d}}){
      print "\n";
      print join("\t", $d, $r)."\n";
      my $data = &process_run($r,$d,$environment);
      push (@bams, $data->{'bam'});
      foreach my $datum (keys(%{$data})){
        print join("\t", $datum, $data->{$datum})."\n";
      }
    }


    my $vcf = $environment->{'output_directory'}."/".$d.".merged.vcf";
    unless (-e $vcf && -f _ && -r _){
      my $mergedrunsbam = $environment->{'output_directory'}."/".$d.".merged.bam";
      unless (-e $mergedrunsbam && -f _ && -r _){
        ## merge all fixed bams
        $cmd = $environment->{'samtools'}." merge ".$mergedrunsbam." ".join(" ", @bams);
        print STDOUT "$cmd\n";
        system($cmd);
      }
      ## run freebayes to generate a per-donor vcf
      $cmd = $environment->{'freebayes'}." -f ".$environment->{'genome_fasta'}." -r NC_000023.11:1-2781479 -@ ".$environment->{'db_snp_vcf'}." -p 1 -0 -b ".$mergedrunsbam." | vcfsnps | vcfbiallelic | vcffilter -f 'QUAL > 20' -f 'TYPE = snp' -f 'AC > 1' -f 'AN > AC + 1' >".$vcf;
      print STDOUT "$cmd\n";
      system($cmd);
    }
    ## get a list of sample RGs (donor, run, cell barcode) that are euploid singlets
    ## call haplotypes on merged vcf with sample list (and banned loci from previous analyses)
  }
}


sub find_executables () {
  my $e;
  $e->{'fastq_dump'};
  $e->{'TagBamWithReadSequenceExtended'} = which('TagBamWithReadSequenceExtended');
  $e->{'samtools'} = which('samtools');
  $e->{'cutadapt'} = which('cutadapt');
  $e->{'bowtie2'} = which('bowtie2');
  $e->{'picard'} = which('picard');
  $e->{'gzip'} = which('gzip');
  $e->{'SpermSeqMarkDuplicates'} = which('SpermSeqMarkDuplicates');
  $e->{'freebayes'} = which('freebayes');
  $e->{'fix_bam_tags'} = which('fix_bam_tags.pl');
  $e->{'haplotypes_from_sperm'} = which('haplotypes_from_sperm.pl');

  return $e;
}


sub process_run ($$$) {
  my ($run, $donor, $environment) = @_;

  unless (defined($environment->{'fastq_directory'})){
    $environment->{'fastq_directory'} = $environment->{'output_directory'};
  }

  #files:
  my $srafile = $environment->{'sra_directory'}."/".$run.".sra";
  my $sraonefqgzfile = $environment->{'fastq_directory'}."/".$run."_1.fastq.gz";
  my $sratwofqgzfile = $environment->{'fastq_directory'}."/".$run."_2.fastq.gz";
  my $unalignedbam = $environment->{'output_directory'}."/".$run.".unal.bam";
  my $cbbam = $environment->{'output_directory'}."/".$run.".unal.cb.bam";
  my $umibam = $environment->{'output_directory'}."/".$run.".unal.cb.umi.bam";
  my $umifqgzfile = $environment->{'output_directory'}."/".$run.".unal.cb.umi.fq.gz";
  my $cutadaptfqgzfile = $environment->{'output_directory'}."/".$run.".unal.cb.umi.cutadapt.fq.gz";
  my $alignedbam = $environment->{'output_directory'}."/".$run.".aligned.cutadapt.bam";
  my $mergebam = $environment->{'output_directory'}."/".$run.".merged.bam";
  my $mkdupbam = $environment->{'output_directory'}."/".$run.".mkdup.bam";
  my $mkdupstats = $environment->{'output_directory'}."/".$run.".mkdup.stats";
  my $fixbam = $environment->{'output_directory'}."/".$run.".fix.bam";
  my $fixtsv = $environment->{'output_directory'}."/".$run.".fix.tsv";
  my $tmpdir = $environment->{'output_directory'}."/tmp";
  my $cmd;

  # processing steps for each run:

  unless (-e $fixbam && -f _ && -r _ && -e $fixtsv && -f _ && -r _){
    unless (-e $mkdupbam && -f _ && -r _ && -e $mkdupstats && -f _ && -r _){
      unless (-e $mergebam && -f _ && -r _){
        unless (-e $alignedbam && -f _ && -r _){
          #unless (-e $cutadaptfqgzfile && -f _ && -r _){
            #unless (-e $umifqgzfile && -f _ && -r _){
              unless (-e $umibam && -f _ && -r _){
                unless (-e $cbbam && -f _ && -r _){
                  unless (-e $unalignedbam && -f _ && -r _){
                    unless (-e $unalignedbam && -f _ && -r _){
                      unless (-e $sraonefqgzfile && -f _ && -r _ && -e $sratwofqgzfile && -f _ && -r _){
                        ## convert SRA to paired fastq.gz
                        $cmd = $environment->{'fastq_dump'}." --gzip  --split-files --outdir ".$environment->{fastq_directory};
                        print STDOUT "$cmd\n";
                        system($cmd);
                      }
                      ## convert paired fastq.gz to unaligned BAM
                      $cmd = $environment->{'picard'}." FastqToSam F1=".$sraonefqgzfile." F2=".$sratwofqgzfile." O=".$unalignedbam." SM=".$donor." RG=".$donor;
                      print STDOUT "$cmd\n";
                      system($cmd);
                    }
                    ## convert paired fastq.gz to unaligned BAM
                    $cmd = $environment->{'picard'}." FastqToSam F1=".$sraonefqgzfile." F2=".$sratwofqgzfile." O=".$unalignedbam." SM=".$donor." RG=".$donor;
                    print STDOUT "$cmd\n";
                    system($cmd);
                  }
                  ## unaligned BAM to new unaligned BAM with cell barcode tag
                  $cmd = $environment->{'TagBamWithReadSequenceExtended'}." BASE_RANGE=1-14 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=true TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 I=".$unalignedbam." O=".$cbbam;
                  print STDOUT "$cmd\n";
                  system($cmd);
                  if (-e $cbbam && -f _ && -r _){
                    unlink $unalignedbam or warn "Could not unlink $unalignedbam: $!";
                  }
                }
                ## unaligned BAM with cell barcode tag to new unaligned BAM with cell barcode tag and UMI tag
                $cmd = $environment->{'TagBamWithReadSequenceExtended'}." BASE_RANGE=1-10 BASE_QUALITY=10 HARD_CLIP_BASES=true BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 I=".$cbbam." O=".$umibam;
                print STDOUT "$cmd\n";
                system($cmd);
                if (-e $umibam && -f _ && -r _){
                  unlink $cbbam or warn "Could not unlink $cbbam: $!";
                }
              }
              ### unaligned BAM with cell barcode tag and UMI tag to single-ended fastq.gz
              #$cmd = $environment->{'samtools'}." fastq -@ 40 -c 9 ".$umibam." | ".$environment->{'gzip'}." -9c - > ".$umifqgzfile;
              #print STDOUT "$cmd\n";
              #system($cmd);
            #}
            ### cutadapt single-ended fastq.gz (min length 20)
            #$cmd = $environment->{'cutadapt'}." -j 40 -a AGATCGGAAGAGC -m 20 -o ".$cutadaptfqgzfile." ".$umifqgzfile;
            #print STDOUT "$cmd\n";
            #system($cmd);

            ### unaligned BAM with cell barcode tag and UMI tag to single-ended fastq.gz with adapters removed
            #$cmd = $environment->{'samtools'}." fastq -@ 40 -c 9 ".$umibam." | ".$environment->{'cutadapt'}." -j 40 -a AGATCGGAAGAGC -m 20 -o ".$cutadaptfqgzfile." - ";
            #print STDOUT "$cmd\n";
            #system($cmd);
          #}
          ### align single-ended fastq.gz
          #$cmd = $environment->{'bowtie2'}." -p 8 --very-sensitive-local -x ".$environment->{'bowtie2_index'}." -U ".$cutadaptfqgzfile." | ".$environment->{'samtools'}." view -b - | ".$environment->{'samtools'}." sort -@ 30 -l 9 -o ".$alignedbam." -";
          #print STDOUT "$cmd\n";
          #system($cmd);

          ## unaligned bam with cb and umi tags to fastq then trim with cutadapt align with bowtie and sort aligned bam
          $cmd = $environment->{'samtools'}." fastq -c 0 ".$umibam." | ".$environment->{'cutadapt'}." -a AGATCGGAAGAGC -m 20 - | ".$environment->{'bowtie2'}." -p 36 --very-sensitive-local -x ".$environment->{'bowtie2_index'}." -U - | ".$environment->{'samtools'}." view -b - | ".$environment->{'samtools'}." sort -l 9 -o ".$alignedbam." -" ;
          print STDOUT "$cmd\n";
          system($cmd);
        }
        ## merge aligned and unaligned bam to recover barcodes and umis
        $cmd = $environment->{'picard'}."  MergeBamAlignment ALIGNED=".$alignedbam." UNMAPPED=".$umibam." O=".$mergebam." R=".$environment->{'genome_fasta'}." INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false TMP_DIR=".$tmpdir;
        print STDOUT "$cmd\n";
        system($cmd);
        if (-e $mergebam && -f _ && -r _){
          unlink $umibam or warn "Could not unlink $umibam: $!";
          unlink $alignedbam or warn "Could not unlink $alignedbam: $!";
        }
      }
      ## mark duplicates by UMI
      # these BAMs are too big to sort in /tmp/
      # option one: use the temp directory
      # $cmd = $environment->{'SpermSeqMarkDuplicates'}." I=".$mergebam." O=".$mkdupbam." CELL_BARCODE_TAG=XC MOLECULAR_BARCODE_TAG=XM NUM_BARCODES=20000 CREATE_INDEX=true OUTPUT_STATS=".$mkdupstats." TMP_DIR=".$tmpdir;
      # option two: use the undocumented option to avoid sorting again!
      # $cmd = $environment->{'SpermSeqMarkDuplicates'}." I=".$mergebam." O=".$mkdupbam." CELL_BARCODE_TAG=XC MOLECULAR_BARCODE_TAG=XM NUM_BARCODES=20000 CREATE_INDEX=true OUTPUT_STATS=".$mkdupstats." SO=unsorted";
      # option three: do both because the SO option still needs a temp directory for some reason
      $cmd = $environment->{'SpermSeqMarkDuplicates'}." I=".$mergebam." O=".$mkdupbam." CELL_BARCODE_TAG=XC MOLECULAR_BARCODE_TAG=XM NUM_BARCODES=20000 CREATE_INDEX=true OUTPUT_STATS=".$mkdupstats." SO=unsorted TMP_DIR=".$tmpdir;
      print STDOUT "$cmd\n";
      system($cmd);
      if (-e $mkdupbam && -f _ && -r _ && -e $mkdupstats && -f _ && -r _){
        unlink $mergebam or warn "Could not unlink $mergebam: $!";
      }
    }
    ## fix bam tags; drop nuisance barcodes; drop low abundance barcodes get new BAM restricted to PAR and a TSV of sex chromosome genotype information
    $cmd = $environment->{'fix_bam_tags'}." -i ".$mkdupbam." -o ".$fixbam." -t ".$fixtsv." -m 14000 -b ".$environment->{'nuisance_barcodes'}." -d ".$donor." -r ".$run;
    print STDOUT "$cmd\n";
    system($cmd);
    if (-e $fixbam && -f _ && -r _ && -e $fixtsv && -f _ && -r _){
      unlink $mkdupbam or warn "Could not unlink $mkdupbam: $!";
    }
  }

  my $run_data;
  $run_data->{'bam'} = $fixbam;
  $run_data->{'tsv'} = $fixtsv;
  $run_data->{'stats'} = $mkdupstats;

  return $run_data;
}


sub usage ($) {
  my $msg = shift(@_);
  my $usage = "$0: $msg\n";
  $usage .= "Usage: $0 -c SraRunTable.csv -s /directory/with/sra_files -b /path/to/bowtie2/index -g genome.fasta -n nuisance.txt -v dbsnp.vcf -o /directory/for/output
\n";
  $usage .= "Try '$0 -h' for more information\n";
  return $usage;
}

sub help (){
  return qq#
$0 (Version: $VERSION)
Usage: $0 -c SraRunTable.csv -s /directory/with/sra_files -b /path/to/bowtie2/index -g genome.fasta -n nuisance.txt -v dbsnp.vcf -o /directory/for/output

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
#;
}
