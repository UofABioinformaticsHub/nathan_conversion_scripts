#!/usr/bin/env perl
# Author: Nathan S. Watson-Haigh (Bioinformatics Hub, University of Adelaide)
# Usage:
#   perl iupac2ab.pl <start_conversion_column> <parent_A_column> <parent_B_column>
use strict;
use warnings;
use Data::Dumper;

my @bases = ('A', 'C', 'G', 'T');

my $first_genotype_col = shift @ARGV;
my $parentA_col = shift @ARGV;
my $parentB_col = shift @ARGV;

my %iupac_ambig = (
  'AG' => 'R',
  'CT' => 'Y',
  'CG' => 'S',
  'AT' => 'W',
  'GT' => 'K',
  'AC' => 'M'
);

while(<>) {
  if ($. == 1) {
    # print the first line (header)
    print;
    next;
  }
  next if /^#/;
  chomp;
  my @columns = split "\t";

  my %tp_bases;
  $tp_bases{$_}++ for @columns[($first_genotype_col-1) .. $#columns];

  if ($columns[$parentA_col-1] =~ /[^ACTG]/) {
    print STDERR "Parent A is not homozygous for the following TagPair: $columns[0] $columns[$parentA_col-1] $columns[$parentB_col-1]\n";
    next;
  } elsif ($columns[$parentB_col-1] =~ /[^ACTG]/) {
    print STDERR "Parent B is not homozygous for the following TagPair: $columns[0] $columns[$parentA_col-1] $columns[$parentB_col-1]\n";
    next;
  } elsif ($columns[$parentA_col-1] eq $columns[$parentB_col-1]) {
    print STDERR "Parent genotypes are the same for the following TagPair: $columns[0] $columns[$parentA_col-1] $columns[$parentB_col-1]\n";
    next;
  }
  
  my %ab_assignment = (
    $columns[$parentA_col-1] => 'AA',
    $columns[$parentB_col-1] => 'BB'
  );
  $ab_assignment{$iupac_ambig{join('', sort($columns[$parentA_col-1], $columns[$parentB_col-1]))}} = 'AB';
 
  #print Dumper(\%ab_assignment); 

  my $line;
  for (my $i=0; $i<($first_genotype_col-1); $i++) {
    $line .= $columns[$i];
    $line .= "\t" unless $i==$first_genotype_col-2;
  }
  #my $line = "$columns[0]\t$columns[1]";
  my $skip_tp=0;
  foreach ( @columns[($first_genotype_col-1) .. $#columns] ) {
    if (/^(N|Z)$/) {
      $line .= "\t$_$_";
    } elsif (exists $ab_assignment{$_}) {
      $line .= "\t$ab_assignment{$_}";
    } elsif (/[ACTG]/) {
      print STDERR "Homozygous individual has a different allele to the parents: $columns[0] $columns[$parentA_col-1] $columns[$parentB_col-1] - $_\n";
      $skip_tp=1;
      last;
    } else {
      print STDERR "Heterozygous individual but different alleles to the parents: $columns[0] $columns[$parentA_col-1] $columns[$parentB_col-1] - $_\n";
      $skip_tp=1;
      last;
    }
  }
  print "$line\n" unless $skip_tp;
}

