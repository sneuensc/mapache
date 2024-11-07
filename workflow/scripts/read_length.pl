#!/usr/bin/perl

# samtools view $bam |length.pl -o ${bam}_length.txt -all
# new: by default N's are ignored (use -all to include N's) 1.11.2024

use Getopt::Long;
my %opts = ();

GetOptions(\%opts, 'o=s', 'all!');

if(scalar(keys(%opts)) < 1){
   &PrintHelp();
}

sub PrintHelp {
   system "pod2text -c $0 ";
   exit();
}

use warnings;
use strict;

my $l = 0;
my @line =();
my %lengths;
my $i = 0;
my $total_bases = 0;
my $average;
my @chromosomes;

open(OUTPUT, ">$opts{'o'}") || die "Cannot create output destination.\n";

while(<STDIN>){
  chomp($_);
  @line = split("\t", $_);

  # Remove 'N' and 'n' from line[9] only if $all is false
  my $filtered_string = $line[9];
  $filtered_string =~ s/[Nn]//g if !$opts{all};

  # Count characters
  my $l = scalar(split("", $filtered_string));
  $total_bases = $total_bases + $l;

  if(exists $lengths{$l}){
    $lengths{$l} = $lengths{$l} + 1;
  }else{
    $lengths{$l} = 1;
  }
  $i = $i + 1;
}

if($i >0){
  $average = $total_bases/$i;
}else{
  $average = "no reads of this type";
}


# print "$average";


print OUTPUT "n_reads\tlength\n";

foreach $i (sort { $a <=> $b } keys %lengths){
  print OUTPUT "$lengths{$i}\t$i\n";
}

close OUTPUT
        or warn $! ? "Error closing : $!"
                   : "Exit status $? ";
