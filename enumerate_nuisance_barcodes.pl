#!/usr/bin/env perl

my $k = 14;
my $l = $k - 1;
my %contaminant;
my %nuisance;
my @bases = ("A", "C", "G", "T");



open(C, "contaminant_list.txt") || die "can't find contaminant_list.txt\n";



while (<C>){
  if (m/^.*\s+([ACGTN]+)$/){
    my $c = $1;
    while (length($c) > $k){
      my $x = substr($c, 0, $k);
      $contaminant{$x} = 1;
      substr($c, 0, 1) = "";
    }
  }
}
close C;

foreach my $c (keys(%contaminant)){
  print $c."\n";
  foreach my $i (0 .. $l){
    my $p = substr($c, $i, 1);
    foreach my $n (@bases){
      unless ($p eq $n){
        my $q = $c;
        substr($q, $i, 1) = $n;
        $nuisance{$q} = 1;
        $nuisance{reverse($q)} = 1;
        $nuisance{complement($q)} = 1;
        $nuisance{reverse_complement($q)} = 1;
      }
    }
  }
}

foreach my $p (@bases){
  my $gcode = $p x $k;
  $nuisance{$gcode} = 1;
  foreach my $i (0 .. $l){
    my $j = $l - $i;
    foreach my $n (@bases){
      unless ($p eq $n){
        if ($i == 0){
          $gcode = $n . $p x $j;
        }else{
          $gcode = $p x $i . $n . $p x $j;
        }
        $nuisance{$gcode} = 1;
      }
    }
  }
}

foreach my $code (keys(%nuisance)){
  print "$code\n";
}

sub complement {
  my $strand = shift(@_);
  $strand =~ tr/ACGT/TGCA/;
  return $strand;
}

sub reverse_complement {
  my $strand = shift(@_);
  my $r = reverse $strand;
  my $rc = complement($r);
  return $rc;
}
