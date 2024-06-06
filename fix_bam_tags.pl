#!/usr/bin/env perl

use Getopt::Long;
use strict;
use vars qw/$VERSION/;

BEGIN {
  $VERSION = '0.2'
};



&main();

sub main {

  ## command line parameters
  my $in_bam_file;
  my $out_bam_file;
  my $rpkm_file;
  my $block_list_file;
  my $donor;
  my $run;
  my $min;
  my $usage;
  my $help;


  #some hard-coded coordinates (sorry):

  my $chrone = {  name => 'chrone',
                  chr => 'NC_000001.11',
                  beg => 1,
                  end => 248956422};


  my $parone = {  name => 'parone',
                  chr => 'NC_000023.11',
                  beg => 1,
                  end => 2781479};

  my $chrnpx = {  name => 'chrnpx',
                  chr => 'NC_000023.11',
                  beg => 2781480,
                  end => 155701382};


  my $partwo = {  name => 'partwo',
                  chr => 'NC_000023.11',
                  beg => 155701382,
                  end => 156030895};

  my $chrnpy = {  name => 'chrnpy',
                  chr => 'NC_000024.10',
                  beg => 2781480,
                  end => 56887902};

  my @regions = ($chrone, $parone, $chrnpx, $partwo, $chrnpy);

  if (!GetOptions(
      'input=s' => \$in_bam_file,
      'output=s' => \$out_bam_file,
      'tsv=s' => \$rpkm_file,
      'blocklist=s' => \$block_list_file,
      'donor=s' => \$donor,
      'run=s' => \$run,
      'minimum=i' => \$min,
      'help' => \$help,
      'usage' => \$usage)){
    print STDERR &usage("$0 (Version: $VERSION)");
    exit -1;
  }

  if ($help){
    print STDOUT &help();
    exit -1;
  }
  if ($usage){
    print STDOUT &usage();
    exit -1;
  }

  my $samtools_exec = `which samtools`;
  chomp($samtools_exec);

  unless (-e $samtools_exec && -x _){
    print STDERR "can't find working samtools executable: $samtools_exec\n";
    exit -1;
  }

  unless (-e $in_bam_file && -r _){
    print STDERR &usage("Please supply an input bam file with -i");
    exit -1;
  }

  unless (defined($out_bam_file)){
    print STDERR &usage("Please supply an output bam file with -o");
  }

  unless (defined($rpkm_file)){
    print STDERR &usage("Please supply an output tsv file with -t");
    exit -1;
  }

  unless (defined($donor)){
    print STDERR "using generic donor id for new read groups\n";
    $donor = "donor";
  }

  unless (defined($run)){
    print STDERR "using generic run id for new read groups\n";
    $run = "run";
  }

  unless (defined($min)){
    $min = 1;
  }
  print STDERR "filtering reads with cell barcodes that have $min or fewer counts\n";

  my %blocklist;
  if (defined($block_list_file)){
    unless (-e $block_list_file && -r _){
      print STDERR &usage("Couldn't read blocklist from: $block_list_file\n");
      exit -1;
    }
    print STDERR "filtering reads with cell barcodes that match blocklist: $block_list_file\n";
    open (B, $block_list_file) || die "Can't open blocklist: $block_list_file\n";
    while (<B>){
      chomp;
      $blocklist{$_} = 1;
    }
    close B;
  }

  #my $in_cmd = "$samtools_exec view -h $in_bam_file | head -n 10000";
  my $in_cmd = "$samtools_exec view -h -F 1024 $in_bam_file";
  my $out_cmd = "$samtools_exec view -h -b - | $samtools_exec sort -l 9 -o $out_bam_file -";

  my %count_rgs;

  open (IN, "$in_cmd |") || die "can't open bam: $in_bam_file with samtools: $samtools_exec\n";
  while (<IN>){
    my $rg; #read_group
    my $xc; #sperm_seq cell_barcode
    if (m/RG:Z:(\S+)/){
      $rg = $1;
    }
    if (m/XC:Z:(\S+)/){
      $xc = $1;
    }
    if ($rg && $xc){
      #print $xc."\t".$blocklist{$xc}."\n";
      unless ($blocklist{$xc}){
        my $new_rg = join("_", $donor, $run, $xc);
        $count_rgs{$new_rg}++;
      }
    }
  }
  close IN;

  my $rg_header = 0;

  my %counts;

  open (OUT, "| $out_cmd") || die "can't write to: $out_bam_file with samtools: $samtools_exec\n";
  open (IN, "$in_cmd |") || die "can't open bam: $in_bam_file with samtools: $samtools_exec\n";
  while (<IN>){
    if (m/^\@/){
      #header lines;
      if (m/^\@RG/){
        #read group header lines;
        unless($rg_header){
          $rg_header = 1;
          foreach my $rg (keys(%count_rgs)){
            unless ($count_rgs{$rg} < $min + 1){
              print OUT join("\t", '@RG', "ID:$rg", "SM:$rg")."\n";
            }
          }
        }
      }else{
        print OUT;
      }
    }else{
      my @sam = split(/\t/, $_);
      my $chr = $sam[2];
      my $pos = $sam[3];
      my $rg; #read_group
      my $xm; #sperm_seq umi
      my $xc; #sperm_seq cell_barcode
      my $cr; #cell_barcode
      my $ox; #umi
      my $mi; #umi

      if (m/RG:Z:(\S+)/){
        $rg = $1;
      }
      if (m/XC:Z:(\S+)/){
        $xc = $1;
        $cr = $xc;
      }
      if (m/XM:Z:(\S+)/){
        $xm = $1;
        $ox = $xm;
        $mi = $xm;
      }

      if ($rg && $xc && $xm){
        my $new_rg = join("_", $donor, $run, $xc);
        #construct new read group name;
        unless ($count_rgs{$new_rg} < $min + 1){
          #skip reads from readgroups with too few reads
          foreach my $r (@regions){
            if ($chr eq $r->{chr}){
              if ($pos >= $r->{beg} && $pos <= $r->{end}){
                #count reads mapping to each region
                $counts{$new_rg}{$r->{name}}++;
                if ($r->{name} eq 'parone'){
                  #only print out PAR reads
                  s/RG:Z:\S+/RG:Z:$new_rg/;
                  s/XC:Z:\S+/XC:Z:$xc\tCR:Z:$cr/;
                  s/XM:Z:\S+/XM:Z:$xm\tOX:Z:$ox\tMI:Z:$mi/;
                  print OUT;
                }
              }
            }
          }
        }
      }
    }
  }
  close IN;
  close OUT;

  open (RPKM, ">$rpkm_file") || die "can't write to TSV: $rpkm_file\n";
  my @names;
  my %kb;
  foreach my $r (@regions){
    push(@names, $r->{name});
    $kb{$r->{name}} = ($r->{end} - $r->{beg} + 1)/1000;
  }
  my $header = join("\t", "RG", "X-bearing", "Sex-aneuploid", "Doublet", @names)."\n";
  print RPKM $header;
  foreach my $rg (keys(%counts)){

    my $million = $count_rgs{$rg}/1000000;
    my $xbearing = 0;
    my $aneuploid = 0;
    my $doublet = 0;
    my %rpkm;
    my @rpkms;

    foreach my $n (@names){
      $rpkm{$n} = ($counts{$rg}{$n}/$kb{$n})/$million;
      push(@rpkms, $rpkm{$n});
    }

    if ($rpkm{'chrnpx'} > $rpkm{'chrnpy'}){
      $xbearing = 1;
    }else{
      $xbearing = 0;
    }

    if ($rpkm{'chrnpx'} > 2 * $rpkm{'chrone'}){
      # supernumerary X chromosome
      # X has more than twice as many reads as chr1
      $aneuploid = 1;
    }elsif ($rpkm{'chrnpy'} > 2 * $rpkm{'chrone'}){
      # supernumerary Y chromosome
      # Y has more than twice as many reads as chr1
      $aneuploid = 1;
    }elsif ($rpkm{'parone'} > 2 * $rpkm{'chrone'}){
      # supernumerary X or Y chromosome; PAR duplications
      # PAR1 has more than twice as many reads as chr1
      # PAR2 is too short and noisy so we will ignore it
      $aneuploid = 1;
    }elsif (2 * $rpkm{'chrnpx'} < $rpkm{'chrone'} && 2* $rpkm{'chrnpy'} < $rpkm{'chrone'}){
      # missing X and Y chromosome
      # X has less than half as many reads as chr1
      # Y has less than half as many reads as chr1
      # PAR1 has less than half as many reads as chr1
      # PAR2 is too short and noisy so we will ignore it
      $aneuploid = 1;
    }else{
      $aneuploid = 0;
    }


    if (
        $rpkm{'chrnpx'} > 2 * $rpkm{'chrnpy'} &&
        $rpkm{'chrnpx'} < 2 * $rpkm{'chrone'} &&
        2 * $rpkm{'chrnpx'} > $rpkm{'chrone'}){
      # euploid X bearing singlet
      # X has more than twice as many reads as Y
      # X has less than twice as many reads as chr1
      # X has more than half as many reads as chr1
      $xbearing = 1;
      $doublet = 0;
    }elsif (
        $rpkm{'chrnpy'} > 2 * $rpkm{'chrnpx'} &&
        $rpkm{'chrnpy'} < 2 * $rpkm{'chrone'} &&
        2 * $rpkm{'chrnpy'} > $rpkm{'chrone'}){
      # euploid Y bearing singlet
      # Y has more than twice as many reads as Y
      # Y has less than twice as many reads as chr1
      # Y has more than half as many reads as chr1
      $xbearing = 0;
      $doublet = 0;
    }else{
      $doublet = 1;
    }

    print RPKM join("\t", $rg, $xbearing, $aneuploid, $doublet, @rpkms)."\n";

  }
  close RPKM;
}

sub usage ($) {
  my $msg = shift(@_);
  my $usage = "$0: $msg\n";
  $usage .= "Usage: $0 -i input.bam -o output.bam -d donor_id -r run_id\n";
  $usage .= "Try '$0 -h' for more information\n";
  return $usage;
}

sub help (){
  return qq#
$0 (Version: $VERSION)
Usage: $0 -i input.bam -o output.bam -t rpkm.tsv

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

#;
}
