#!/usr/bin/env perl

if ($#ARGV < 0) {
  print "usage: $0 infile [-c]\n";
  print "-c will print the contig name\n";
  exit(0);
}

$infile = shift @ARGV;
$printContig = 0;
if ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-c") {
    $printContig = 1;
  }
}
open(IN, "$infile") or die "cannot open $infile\n";

$contigSize = -1;
$contig = "";
$first = 1;
while(<IN>) {
  $line = $_;
  chomp($line);
  if (/^>/) {
    if ($first == 0) {
				if ($contigSize > 0 ) {
						if ($printContig == 1) {
								print "$contigSize $contig\n";
						} 
						else {
								print "$contigSize\n";
						}
				}
    }
    $first = 0;
    $contigSize = 0;
		$line =~ />(\S+)/;
    $contig = $1;
  }
  else {
    $l = length($line);
    $contigSize += length($line);
  }
}

if ($contigSize > 0 ) {
  if ($printContig == 1) {
    print "$contig $contigSize\n";
  } 
  else {
    print "$contigSize\n";
  }
}
