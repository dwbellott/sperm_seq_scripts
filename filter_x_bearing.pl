#!/usr/bin/env perl

while (<STDIN>){
  my $x = $_;
  chomp($x);
  my @line = split(/\t/, $x);
  my ($chrone, $parone, $chrnpx, $partwo, $chrnpy) = ($line[4], $line[5], $line[6], $line[7], $line[8]);
  if ($chrnpx < 2* $chrone && $chrnpx > 0.5* $chrone){
    if ($chrnpy < 0.25 * $chrone && $chrnpy < 0.5 * $chrnpx){ 
      print;
    }
  }  
}
