#!/usr/bin/env perl
use Getopt::Long;
use strict;
use GD::SVG;
use vars qw/$VERSION/;

BEGIN {
  $VERSION = '0.5'
};


## command line parameters
my $vcf_file;
my $chromosome_name;
my $karyotype_file;
my $sample_file


&main();


sub main {

  #setup command-line options
  my (
    $vcf_file,
    $chromosome_name,
    $npx_position,
    $karyotype_file,
    $sample_file,
    $xnco_file,
    $ynco_file,
    $ban_file,
    $image_output,
    $co_location_output,
    $bed_graph_output,
    $help,
    $usage);

  if (!GetOptions(
      'vcf=s' => \$vcf_file,
      'chromosome=s' => \$chromosome_name,
      'position=i' => \$npx_position,
      'karyotype=s' => \$karyotype_file,
      'samples=s' => \$sample_file,
      'xnco=s'  => \$xnco_file,
      'ynco=s'  => \$ynco_file,
      'banned=s'  => \$ban_file,
      'image=s' => \$image_output,
      'location=s' => \$co_location_output,
      'graph=s' => \$bed_graph_output,
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

  unless(defined($vcf_file)){
    print STDERR &usage("please supply a VCF file with -v");
    exit -1;
  }

  unless(defined($chromosome_name)){
    print STDERR &usage("please supply a chromosome name with -c");
    exit -1;
  }

  unless(defined($npx_position)){
    print STDERR &usage("please supply the start position of the NPX/NPY with -p");
    exit -1;
  }

  unless(defined($karyotype_file)){
    print STDERR &usage("please supply a karyotype file with -k");
    exit -1;
  }

  unless(defined($sample_file)){
    print STDERR &usage("please supply a list of samples with -s");
    exit -1;
  }

  unless(defined($image_output)){
    $image_output = "output.image.svg";
  }

  unless(defined($co_location_output)){
    $co_location_output = "output.locations.tsv";
  }

  unless(defined($bed_graph_output)){
    $bed_graph_output = "output.recombination.bedGraph";
  }

  &run_analysis($vcf_file, $chromosome_name, $karyotype_file, $sample_file, $image_output, $co_location_output, $npx_position, $bed_graph_output);
  exit 0;
}

sub usage ($) {
  my $msg = shift(@_);
  my $usage = "$0: $msg\n";
  $usage .= "Usage: $0 -v variants.vcf -k karyotype.tsv -s samples.txt -c chomosome_name\n";
  $usage .= "Try '$0 -h' for more information\n";
  return $usage;
}

sub help (){
  return qq#
$0 (Version: $VERSION)
Usage: $0 -v input.vcf -t calls.tsv -s samples.txt -c chrX -p 2781480

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

#;
}

## subroutines

sub run_analysis{
  my ($vcf_file, $chromosome_name, $karyotype_file, $sample_file, $image_output, $co_location_output, $npx_position, $bed_graph_output) = @_;

  my $five_centimorgan_distance = (5 / 50) * $npx_position;
  # no filter to start
  my $filter;
  $filter->{filtering} = 0;
  #get variants for our chromosome;
  print STDERR "begin first pass\n";
  my $variants = &parse_vcf($vcf_file, $chromosome_name, $filter);
  #get sex chromosome genotype;
  my $x_chr;
  if (defined($karyotype_file)){
    $x_chr = &parse_karyotypes($karyotype_file);
    $x_chr->{chromosome} = $chromosome_name;
    $x_chr->{position} = $npx_position;
    push (@{$variants}, $x_chr);
  }
  # load samples
  my $samples = &load_samples($sample_file);
  # perform initial haplotyping on all variants
  my $haplotype = &reconstruct_haplotypes($variants, $samples);
  # use haplotyping to identify crossovers
  my $crossovers = &identify_crossovers($haplotype, $variants, $samples);
  # use crossovers to identify unreliable markers
  my $unreliable = &identify_unreliable_markers($crossovers, $variants, $samples, $five_centimorgan_distance);
  # unreliable markers go to: filter->banned;
  my $filter->{banned} = &ban_unreliable_markers($unreliable);
  # crossovers involving unreliable markers are filtered
  my $filtered_xo = &filter_unreliable_crossovers($crossovers, $unreliable);
  # filtered crossovers
  # identify NCO samples
  my $nco_samples = &identify_nco_samples($crossovers, $variants, $samples);
  # filter variants based on NCO samples
  $filter->{x} = $nco_samples->{x};
  $filter->{y} = $nco_samples->{y};
  $filter->{filtering} = 1;
  #get variants for our chromosome;
  print STDERR "begin second pass\n";
  my $variants2 = &parse_vcf($vcf_file, $chromosome_name, $filter);
  #get sex chromosome genotype;
  if (defined($x_chr)){
    $x_chr->{chromosome} = $chromosome_name;
    $x_chr->{position} = $npx_position;
    push (@{$variants2}, $x_chr);
  }
  # perform second round of haplotyping on filtered variants
  my $haplotype2 = &reconstruct_haplotypes($variants2, $samples);
  # use 2nd haplotyping to identify crossovers
  my $crossovers2 = &identify_crossovers($haplotype2, $variants2, $samples);
  # sort samples
  my $filtered_xo2 = &filter_tight_double_recombinants($crossovers2,$five_centimorgan_distance);
  my $samples_sorted = &sort_samples($filtered_xo2, $variants2, $samples);
  # produce an image?
  if ($image_output){
    &print_svg($image_output, $haplotype2, $variants2, $samples_sorted, $filtered_xo2);
  }
  if ($co_location_output){
    &print_cos($co_location_output, $filtered_xo2);
  }
  if ($bed_graph_output){
    &print_bed_graph($bed_graph_output, $filtered_xo2, $samples_sorted, $variants2);
  }
  return 1;
}

sub filter_tight_double_recombinants {
  my ($crossovers, $five_centimorgan_distance) = @_;

  my $return_xo;

  my %tight;
  my $last_xo;

  foreach my $xo (@{$crossovers}){
    if (  defined($last_xo) &&
          $xo->{sample} eq $last_xo->{sample} &&
          $xo->{previous}->{chromosome} eq $last_xo->{next}->{chromosome}){

      # double crossover situation
      # same sample
      my $s = $xo->{sample};
      # same chromosome
      my $c = $xo->{previous}->{chromosome};

      my $window = $xo->{next}->{position} - $last_xo->{previous}->{position};

      #print STDERR join("\t", $xo->{sample}, $c, $last_xo->{previous}->{position}, $last_xo->{next}->{position}, $xo->{previous}->{position}, $xo->{next}->{position}, $window, $five_centimorgan_distance);

      if ($last_xo->{next}->{position} == $xo->{previous}->{position}){
        # if the dco relies on a single marker
        # probably a genotyping error
        # invalidate both flanking crossovers
        $tight{$s}{$c}{$last_xo->{previous}->{position}}{$last_xo->{next}->{position}} = 1;
        $tight{$s}{$c}{$xo->{previous}->{position}}{$xo->{next}->{position}} = 1;
        #print STDERR "\t**genotyping error\n"
      }elsif ($window < $five_centimorgan_distance){
        # if the dco window is narrower than 5cM
        # probably gene conversion
        # invalidate both flanking crossovers
        $tight{$s}{$c}{$last_xo->{previous}->{position}}{$last_xo->{next}->{position}} = 1;
        $tight{$s}{$c}{$xo->{previous}->{position}}{$xo->{next}->{position}} = 1;
        #print STDERR "\t**gene conversion\n"
      }else{
        #print STDERR "\n";
      }
    }
    $last_xo = $xo;
  }

  foreach my $xo (@{$crossovers}){
    if ($tight{$xo->{sample}}{$xo->{previous}->{chromosome}}{$xo->{previous}->{position}}{$xo->{next}->{position}}){
      #print STDERR join("\t", $xo->{sample}, $xo->{previous}->{chromosome}, $xo->{previous}->{position}, $xo->{next}->{position}, $five_centimorgan_distance)."\n";
    }else{
      push(@{$return_xo}, $xo);
    }
  }
  return $return_xo;
}

sub load_samples{
  my $sample_file = shift(@_);
  my $samples;
  open (S, $sample_file) || die "can't read sample names from $sample_file\n";
  while (<S>){
    chomp;
    push(@{$samples}, $_);
  }
  return $samples;
}

sub parse_karyotypes{
  my $karyotype = shift(@_);
  my $x_chr;
  open (K, $karyotype) || die "can't read calls in $karyotype\n";
  while (<K>){
    chomp;
    my @tsv = split(/\t/, $_);
    my $sm = $tsv[0];
    my $xb = $tsv[1];

    #print "$sm\t$xb\n";
    $x_chr->{$sm} = $xb;
  }
  close K;
  return $x_chr;
}

sub parse_vcf {
  my ($vcf_file, $chromosome, $filter) = @_;
  my $variants;

  my %column_read_group;
  my @read_groups;

  my $last_position;

  open(VCF, $vcf_file) || die "can't read from $vcf_file\n";
  while(<VCF>){
    if (m/^\#/){
      if (m/^#CHROM/){
        chomp;
        my @vcf = split("\t", $_);
        foreach my $i (9 .. $#vcf){
          $column_read_group{$i} = $vcf[$i];
          push(@read_groups, $vcf[$i]);
        }
      }
    }else{
      my @vcf = split("\t", $_);
      my $seven = &parse_seven($vcf[7]);
      if (  $vcf[0] eq $chromosome &&
            !defined($filter->{banned}->{$vcf[1]})){
        #on the correct chromosome and not banned
        if (  $vcf[5] >= 30 &&
              $seven->{"TYPE"} eq "snp" &&
              $seven->{"AC"} > 1 &&
              $seven->{"AN"} > $seven->{"AC"} + 1 &&
              $seven->{"AF"} > 0.25 &&
              $seven->{"AF"} < 0.75){
          #Filter variant calls:
          # high quality calls
          # snp markers only
          # at least two each ref and alt alleles
          # observed allele frequency between 0.25 and 0.75

          my $site;

          if ($filter->{filtering} == 1){
            # avoid undefs
            $site->{x}->{0} = 0;
            $site->{x}->{1} = 0;
            $site->{y}->{0} = 0;
            $site->{y}->{1} = 0;
          }

          $site->{chromosome} = $vcf[0];
          $site->{position} = $vcf[1];
          $site->{quality} = $vcf[5];
          $site->{depth} = $seven->{"DP"};
          foreach my $i (9 .. $#vcf){
            my @fields = split(":", $vcf[$i]);
            my ($gt,$dp) = ($fields[0],$fields[1]);
            $site->{$column_read_group{$i}} = $gt;
            if ($filter->{x}->{$column_read_group{$i}}){
              $site->{x}->{$gt}++;
            }elsif ($filter->{y}->{$column_read_group{$i}}){
              $site->{y}->{$gt}++;
            }
          }

          if ($filter->{filtering} == 1){
            # if we're filtering
            # take only consistent markers

            my $consistency = 0;

            #case 1: we have no allele data for x and y samples:

            if (  $site->{x}->{0} == 0 &&
                  $site->{x}->{1} == 0 &&
                  $site->{y}->{0} == 0 &&
                  $site->{y}->{1} == 0){
              #case 0: we have no allele data for x and y samples:
              # we don't want this marker
              $consistency = 0;
            }elsif (  $site->{x}->{0} == 0 &&
                      $site->{x}->{1} == 0){
              #case 1: we have no allele data for x samples:
              # we could take these if the y samples are consistent:
              # $consistency = $site->{y}->{0} / ($site->{y}->{0} + $site->{y}->{1});
              # but we won't right now:
              $consistency = 0;
            }elsif (  $site->{y}->{0} == 0 &&
                      $site->{y}->{1} == 0){
              #case 2: we have no allele data for y samples:
              # we could take these if the x samples are consistent:
              # $consistency = $site->{x}->{0} / ($site->{x}->{0} + $site->{x}->{1});
              # but we won't right now:
              $consistency = 0;
            }else{
              #case 2: we have allele data for x and y samples:
              my $xreffraction = $site->{x}->{0} / ($site->{x}->{0} + $site->{x}->{1});
              my $yreffraction = $site->{y}->{0} / ($site->{y}->{0} + $site->{y}->{1});
              # these should be inversely related
              $consistency = abs($xreffraction - $yreffraction);
            }

            if ($consistency > 0.9){
              push (@{$variants}, $site);
              $last_position = $vcf[1];
            }

          }elsif ($vcf[1] > $last_position + 1){
            push (@{$variants}, $site);
            $last_position = $vcf[1];
          }

        }
      }
    }
  }
  close VCF;
  return $variants
}

sub parse_seven {
  my $field = shift(@_);
  my $result;
  my @seven = split(/\;/, $field);
  foreach my $s (@seven){
    my ($tag,$value) = split(/\=/, $s);
    $result->{$tag} = $value;
  }
  return $result;
}

sub reconstruct_haplotypes {
  my ($variants, $samples) = @_;
  my $sites = $#{$variants};
  my $quarter = int($sites/4)+1;
  my $haplotype;
  $haplotype->[$sites]->{0} = 0;
  $haplotype->[$sites]->{1} = 1;
  foreach my $i (reverse(0 .. $sites -1)){
    if ($i < 2800){
      #die;
    }
    my ($n, $m, $p, $q) = (0, 0, 0, 0);
    foreach my $j (reverse($i + 1 .. min($sites, $i + $quarter))){
      foreach my $s (@{$samples}){
        if (defined($haplotype->[$j]->{0}) && $variants->[$i]->{$s} ne "." && $variants->[$j]->{$s} ne "."){
          #print STDERR "i: $i\tj: $j\ts: $s\tsi: $variants->[$i]->{$s}\tsj: $variants->[$j]->{$s}\thj0: $haplotype->[$j]->{0}\thj1: $haplotype->[$j]->{1}\n";
          if ($variants->[$i]->{$s} eq "0" && $variants->[$j]->{$s} eq $haplotype->[$j]->{0}){
            $n++;
          }elsif ($variants->[$i]->{$s} eq "0" && $variants->[$j]->{$s} eq $haplotype->[$j]->{1}){
            $m++;
          }elsif ($variants->[$i]->{$s} eq "1" && $variants->[$j]->{$s} eq $haplotype->[$j]->{0}){
            $p++;
          }elsif ($variants->[$i]->{$s} eq "1" && $variants->[$j]->{$s} eq $haplotype->[$j]->{1}){
            $q++;
          }
        }
      }
    }
    #print STDERR "site $i: 0/0 = $n 0/1 = $m 1/0 = $p 1/1 = $q\n";
    if ($n + $q > $m + $p){
      $haplotype->[$i]->{0} = 0;
      $haplotype->[$i]->{1} = 1;
      #print STDERR "site $i: h0 = 0 and h1 = 1\n";
    }elsif ($n + $q < $m + $p){
      $haplotype->[$i]->{0} = 1;
      $haplotype->[$i]->{1} = 0;
      #print STDERR "site $i: h0 = 1 and h0 = 0\n";
    }else{
      #print STDERR "could not phase variant site $i\n";
    }
  }
  return $haplotype;
}

sub identify_crossovers {
  my ($haplotype, $variants, $samples) = @_;
  my $sites = $#{$variants};
  my $crossovers;
  foreach my $s (@{$samples}){
    my $previous;
    foreach my $i (0 .. $sites){
      my $current;
      $current->{chromosome} = $variants->[$i]->{chromosome};
      $current->{position} = $variants->[$i]->{position};
      #print join("\t", $s, $i, $variants->[$i]->{$s}, $haplotype->[$i]->{1}, $haplotype->[$i]->{0})."\n";
      if ($variants->[$i]->{$s} eq "."){
        undef $current->{haplotype};
      }elsif ($variants->[$i]->{$s} == $haplotype->[$i]->{1}){
        $current->{haplotype} = "X";
      }elsif ($variants->[$i]->{$s} == $haplotype->[$i]->{0}){
        $current->{haplotype} = "Y";
      }else{
        undef $current->{haplotype};
      }
      if (defined($current->{haplotype})){
        if (defined($previous->{haplotype})){
          unless ($previous->{haplotype} eq $current->{haplotype}){
            #print join("\t", $s,$previous->{chromosome}, $previous->{position}, $previous->{haplotype}, $current->{chromosome}, $current->{position}, $current->{haplotype})."\n";
            my $xo;
            $xo->{sample} = $s;
            $xo->{previous} = $previous;
            $xo->{next} = $current;
            push(@{$crossovers}, $xo);
          }
        }
        $previous = $current;
      }
    }
  }
  return $crossovers;
}

sub identify_unreliable_markers {
  my ($crossovers, $variants, $samples, $five_centimorgan_distance) = @_;
  my $sites = $#{$variants};
  my $last_xo;
  my %count;
  my %appearances;

  #print STDERR "\tlooking for unreliable markers\n";


  foreach my $i (0..$sites){
    my ($c,$p) = ($variants->[$i]->{chromosome},$variants->[$i]->{position});
    foreach my $s (@{$samples}){
      if ($variants->[$i]->{$s} eq "."){
      }elsif ($variants->[$i]->{$s} eq "1"){
          $appearances{$c}{$p}++;
      }elsif ($variants->[$i]->{$s} eq "0"){
          $appearances{$c}{$p}++;
      }
    }
  }

  foreach my $xo (@{$crossovers}){
    if (  defined($last_xo) &&
          $xo->{sample} eq $last_xo->{sample} &&
          $xo->{previous}->{chromosome} eq $last_xo->{next}->{chromosome}){

      # double crossover situation
      # same sample
      # same chromosome

      my $c = $xo->{previous}->{chromosome};
      my $window = $xo->{next}->{position} - $last_xo->{previous}->{position};
      #my $window = $xo->{previous}->{position} - $last_xo->{next}->{position};

      #print STDERR join("\t", $xo->{sample}, $c, $last_xo->{previous}->{position}, $last_xo->{next}->{position}, $xo->{previous}->{position}, $xo->{next}->{position}, $window, $five_centimorgan_distance)."\n";

      if ($last_xo->{next}->{position} == $xo->{previous}->{position}){
        # if the dco relies on a single marker
        # probably a genotyping error
        # mark the position as unreliable
        my $p = $xo->{previous}->{position};
        #print STDERR "**\t$c\t$p\n";
        $count{$c}{$p}++;
      }elsif ($window < $five_centimorgan_distance){
        # if the dco window is narrower than 5cM
        # probably gene conversion
        # and mark intervening positions as unreliable
        foreach my $p ($last_xo->{next}->{position} .. $xo->{previous}->{position}){
          if ($appearances{$c}{$p}){
            #print STDERR "##\t$c\t$p\n";
            $count{$c}{$p}++;
          }
        }
      }
    }
    $last_xo = $xo;
  }

  my $unreliable;
  foreach my $c (keys(%count)){
    foreach my $p (keys(%{$count{$c}})){
      my $u;
      $u->{chromosome} = $c;
      $u->{position} = $p;
      $u->{count} = $count{$c}{$p};
      $u->{appearances} = $appearances{$c}{$p};
      $u->{ratio} = $u->{count}/$u->{appearances};
      push(@{$unreliable}, $u);
    }
  }
  return $unreliable;
}

sub ban_unreliable_markers {
  my ($unreliable) = @_;
  my $filter;
  foreach my $u (@{$unreliable}){
    if ($u->{ratio} > 0.001){
      #we're aiming for Q=30 which is 1 error in 1000
      $filter->{$u->{position}} = 1;
    }
  }
  return $filter;
}

sub filter_unreliable_crossovers {
  my ($crossovers, $unreliable) = @_;
  my $filtered;
  my %filter;
  foreach my $u (@{$unreliable}){
    $filter{$u->{chromosome}}{$u->{position}} = 1;
  }
  foreach my $xo (@{$crossovers}){
    if ($filter{$xo->{previous}->{chromosome}}{$xo->{previous}->{position}}){
      #do nothing
    }elsif ($filter{$xo->{next}->{chromosome}}{$xo->{next}->{position}}){
      #do nothing
    }else{
      push(@{$filtered}, $xo);
    }
  }
  return $filtered;
}


sub identify_nco_samples {
  my ($crossovers, $variants, $samples) = @_;

  my $nco;

  #check chromosome constitution
  my $x_chr = $variants->[$#{$variants}];

  my %count;

  foreach my $s (@{$samples}){
    $count{$s} = 0;
  }

  foreach my $xo (@{$crossovers}){
    $count{$xo->{sample}}++;
  }

  foreach my $s (@{$samples}){
    if ($count{$s} == 0){
      if ($x_chr->{$s} == 1){
        $nco->{x}->{$s} = 1;
      }else{
        $nco->{y}->{$s} = 1;
      }
    }
  }

  return $nco;

}

sub sort_samples {
  my ($crossovers, $variants, $samples) = @_;
  my %next;
  my %prev;
  my %count;

  foreach my $s (@{$samples}){
    $count{$s} = 0;
  }

  foreach my $xo (@{$crossovers}){
    $count{$xo->{sample}}++;
    $next{$xo->{sample}} = $xo->{next}->{position};
    $prev{$xo->{sample}} = $xo->{previous}->{position};
  }

  my $x_chr = $variants->[$#{$variants}];
  my $sorted;
  @{$sorted} = sort {
    $x_chr->{$b} <=> $x_chr->{$a} ||
    $count{$a} <=> $count{$b} ||
    $next{$a} <=> $next{$b} ||
    $prev{$a} <=> $prev{$b} ||
    $a cmp $b

    } @{$samples};

  return $sorted;
}

sub print_cos {
  my ($co_location_output, $crossovers) = @_;
  open (COS, ">$co_location_output") || die "can't write crossovers to $co_location_output\n";
  foreach my $xo (@{$crossovers}){
    print COS join("\t", $xo->{sample}, $xo->{previous}->{chromosome}, $xo->{previous}->{position}, $xo->{next}->{position})."\n";
  }
  close COS;
  return 1;
}

sub print_bed_graph {
  my ($bed_graph_output, $crossovers, $samples, $variants) = @_;

  my $sample_count = $#{$samples};
  my $individual_rf = (1 / $sample_count) * 100;
  my $last_variant = $#{$variants};

  my $chrom = $variants->[$last_variant]->{chromosome};
  my $first_pos = $variants->[0]->{position};
  my $last_pos = $variants->[$last_variant]->{position};

  my @cm_per_mb = (0) x $last_pos;

  foreach my $xo (@{$crossovers}){
    my $beg = $xo->{previous}->{position};
    my $end = $xo->{next}->{position};
    my $length_in_mb = ($end - $beg) / 10**6;
    my $local_rate = $individual_rf / $length_in_mb;
    foreach my $i ($beg .. $end){
      $cm_per_mb[$i] += $local_rate;
    }
  }

  my $bed;

  open (BEDGRAPH, ">$bed_graph_output") || die "can't write bedGrah to $bed_graph_output\n";

  print BEDGRAPH 'track type=bedGraph name="Local Recombination Rate" description="Local Recombination Rate Calculated from SpermSeq Crossovers across '.$sample_count.' Sperm" visibility=full color=200,100,0 altColor=0,100,200 priority=20'."\n";

  foreach my $i (1 .. $last_pos){
    if (defined($bed->{chromStart})){
      if ($bed->{dataValue} == $cm_per_mb[$i]){
        #extend current line
        $bed->{chromEnd} = $i;
      }else{
        print BEDGRAPH join("\t", $chrom, $bed->{chromStart}, $bed->{chromEnd}, $bed->{dataValue})."\n";
        #start fresh line
        $bed->{chromStart} = $i;
        $bed->{chromEnd} = $i;
        $bed->{dataValue} = $cm_per_mb[$i];
      }
    }else{
      #start fresh line
      $bed->{chromStart} = $i;
      $bed->{chromEnd} = $i;
      $bed->{dataValue} = $cm_per_mb[$i];
    }
  }

}

sub print_svg{
  my ($image_file, $haplotype, $variants, $samples, $crossovers) = @_;

  my $sites = $#{$variants};
  my $base_width = $variants->[$sites]->{position};
  my $chromosome_count = $#{$samples} + 1;

  my $px_margin = 20;
  my $px_paper_width = 2551; #US letter
  my $px_min_tick = 2;
  my $px_per_base = ($px_paper_width-2*$px_margin-$px_min_tick)/$base_width;

  my $tick_width = $px_per_base * 1;
  if ($tick_width < 1){
    $tick_width = $px_min_tick;
  }

  my $chromosome_width = $px_per_base * $base_width;
  my $chromosome_height = 10;

  my $width = $px_paper_width;
  my $height = 3*$px_margin + $chromosome_count*($px_margin + $chromosome_height);

  my $img = GD::SVG::Image->new($width, $height);
  $img->interlaced('true');

  my $white = $img->colorAllocate(255,255,255);
  my $black = $img->colorAllocate(0,0,0);
  my $grey = $img->colorAllocate(187,187,187);
  my $ltgrey = $img->colorAllocate(204,204,204);
  my $dkgrey = $img->colorAllocate(153,153,153);
  my $pargreen = $img->colorAllocate(179,211,190);
  my $x_orange = $img->colorAllocate(222,112,38);
  my $y_purple = $img->colorAllocate(120,83,162);
  my $red = $img->colorAllocate(255, 0, 0);
  my $blue = $img->colorAllocate(0, 0, 255);

  my $y = $px_margin;
  my $x = $px_margin;

  my $p = 0;
  my $ytext = $y;
  my $yline0 = $y + 2*$px_margin;
  my $yline1 = $yline0 + $chromosome_count*($px_margin + $chromosome_height) - $px_margin;


  $y += 2*$px_margin;

  foreach my $s (@{$samples}){

      $img->string(gdSmallFont,$x,$y-($px_margin/2),$s,$black);

      $img->filledRectangle($x,$y,$x+$chromosome_width,$y+$chromosome_height,$pargreen);

      ### depict crossovers:
      foreach my $xo (@{$crossovers}){
        print COS join("\t", $xo->{sample}, $xo->{previous}->{chromosome}, $xo->{previous}->{position}, $xo->{next}->{position})."\n";
        if ($xo->{sample} eq $s){
          my $x0 = $x + $xo->{previous}->{position} * $px_per_base;
          my $x1 = $x + $xo->{next}->{position} * $px_per_base;
          my $y0 = $y;
          my $y1 = $y + $chromosome_height;
          $img->filledRectangle($x0,$y0,$x1,$y1,$dkgrey);
        }
      }

      ### depict markers:
      foreach my $i (0..$sites){

        my $x0 = $x + $variants->[$i]->{position} * $px_per_base;
        my $x1 = $x0 + $tick_width;
        my $y0 = $y;
        my $y1 = $y + $chromosome_height;

        if ($haplotype->[$i]->{1} eq $variants->[$i]->{$s}){
          $img->filledRectangle($x0,$y0,$x1,$y1,$x_orange);
        }elsif ($haplotype->[$i]->{0} eq $variants->[$i]->{$s}){
          $img->filledRectangle($x0,$y0,$x1,$y1,$y_purple);
        }

      }

      $y = $y + $chromosome_height + $px_margin;
  }

  while ($p < $base_width){
    my $xline = $x + $p * $px_per_base;
    if ($p % 1e6 == 0){
      my $string = $p / 1e6;
      $img->string(gdMediumBoldFont,$xline,$ytext,$string,$black)
    }
    $img->line($xline,$yline0,$xline,$yline1,$white);
    $p += 250000;
  }
  $p = $base_width - 500000;
  my $xline = $x + $p * $px_per_base;
  $img->line($xline,$yline0,$xline,$yline1,$black);

  open(SVG, ">$image_file") || die "can't write svg image: $image_file\n";
  print SVG $img->svg;
  close SVG;
}

sub max {
  my ($max, @others) = @_;
  while (@others){
    my $next = shift(@others);
    if ($next > $max){
      $max = $next;
    }
  }
  return $max;
}

sub min {
  my ($min, @others) = @_;
  while (@others){
    my $next = shift(@others);
    if ($next < $min){
      $min = $next;
    }
  }
  return $min;
}
