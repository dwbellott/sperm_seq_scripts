#!/usr/bin/env perl

my @beds = @ARGV;

my $individuals = $#beds;

my %sum;

my $max;

foreach my $bedfile (@beds){
  open (BED, $bedfile) || die "can't read $bedfile\n";
  while (<BED>){
    chomp;
    my ($chr, $beg, $end, $data) = split($_);
    if ($max < $end){
      $max = $end;
    }
    foreach my $i ($beg .. $end){
      $sum{$i} += $data;
    }
  }
}

print 'track type=bedGraph name="Local Recombination Rate"'."\n";

my @cm_per_mb = (0) x $max;

foreach my $x (keys(%sum)){
  $cm_per_mb[$x] = $sum{$x}/$individuals;
}

my $bed;

foreach my $i (1 .. $max){
  if (defined($bed->{chromStart})){
    if ($bed->{dataValue} == $cm_per_mb[$i]){
      #extend current line
      $bed->{chromEnd} = $i;
    }else{
      print join("\t", $chrom, $bed->{chromStart}, $bed->{chromEnd}, $bed->{dataValue})."\n";
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
