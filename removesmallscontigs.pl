#!/usr/bin/perl

use strict;
use warnings;

# Usage: perl removesmalls.pl 1000 input.fasta > output.fasta

my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";  # change temporarely the input record separator from "\n" to ">"
    while(<>) {
        chomp;
        next unless /\w/;  # skip line if no word character ("\w") is found
        s/>$//gs;  # 
        my @chunk = split /\n/;
        my $header = shift @chunk;  # first item in array
        my $seqlen = length join "", @chunk;  # all the other item of the array
        print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";  #restore input record separator
}