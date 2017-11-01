#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;  # debug

#I\O files
my $input_sampleList;  # a text files with all samples listed. One sample per line.
my $input_clstr;  # Output ".clstr" file from CD-HIT-EST
my $output_table;  # Tab-separated output table
my $help;
my $man;

GetOptions(
    'sample_list|s=s' => \$input_sampleList,
    'input_clstr|i=s' => \$input_clstr,
    'output_table|o=s' => \$output_table,
    'help|h' => \$help,
    'man' => \$man) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;


#check if arguments are not empty
if (!defined($input_sampleList) || !defined($input_clstr) || !defined($output_table))
{
    print ("Not all required arguments supplied.\nUsage: perl cdHitClstr2table.pl -s sampleList -c cd-hit.clstr -o outputTable.tsv\n");
    exit 1;
#    die "Usage: perl cdHitClstr2table.pl -s sampleList -c cd-hit.clstr -o outputTable.tsv";
}


##############################
#                            #
#   Parse sample list file   #
#                            #
##############################


# open sample list file in read-only mode
open (my $input_sampleList_fh, "<", $input_sampleList) or die "Cannot open $input_sampleList: $!\n";

my @samples;

# Read file line by line
while (my $line = <$input_sampleList_fh>)
{
    chomp($line);  # remove carriage return at end of line
    next if $line eq "";  # skip empty lines
    push (@samples, $line);
}

close($input_sampleList_fh);


###############################
#                             #
#   Parse cd-hit-est report   #
#                             #
###############################


# Open cd-hit-est cluster report file in read-only mode
open (my $input_clstr_fh, "<", $input_clstr) or die "Cannot open $input_clstr: $!\n";

# Put sequence from list in an array
my %clusters;  # Cluster ID
my $clust;

# Read line by line
while (my $line = <$input_clstr_fh>)
{
    chomp($line);  # remove carriage return at end of line
    next if $line eq "";  # skip empty lines

    # >Cluster 0
    # 0 78657nt, >2016-SEQ-0017_2... *
    # 1 78657nt, >2016-SEQ-0013_2... at +/100.00%

    if ($line =~ /^>Cluster/)
    {
        # Name of the cluster
        $line =~ s/>//;  # remove ">" (substitute by nothing)
        $line =~ s/ /_/;  # substitute space by underscore
        $clust = $line;
    }
    else
    {
        my @fields = split(/\t/, $line);  # There are 2 tab-separated fields
        my @info = split(/ /, $fields[1]);  # split at space

        # Sample ID
        my $id = (split(/_/, $info[1]))[0];
        $id =~ s/\.//g;
        $id =~ s/>//;

        # Length of members
        my $len = $info[0];
        $len =~ s/nt\,//;  # only keep the numerical value

        # Add length to hash
        $clusters{$clust}{$id} = $len;
    }
}

close ($input_clstr_fh);

# print Dumper \%clusters;


##############################
#                            #
#   Write ouput table file   #
#                            #
##############################


# Open output file to write
open (my $output_table_fh, ">", $output_table) or die "Cannot open $output_table: $!";

# Print
# print {$output_table_fh} ("# Constructed from biom file" . "\n");

# Print header with sample names
print {$output_table_fh} ("#OTU ID" . "\t" . join("\t", sort(@samples)) . "\n");

# Loop through hash to print output table
foreach my $cluster (sort(keys(%clusters)))
{
    print {$output_table_fh} ($cluster);

    foreach my $sample (sort(@samples))  # keep the same order than the header
    {
        if ($clusters{$cluster}{$sample})
        {
            print {$output_table_fh} ("\t" . $clusters{$cluster}{$sample});
        }
        else
        {
            print {$output_table_fh} ("\t" . "0");
        }
    }
    print {$output_table_fh} ("\n");
}

close($output_table_fh);

__END__
=head1 NAME

cdHitClstr2table.pl - Convert CD-HIT-EST ".clstr" output file for QIIME.

=head1 SYNOPSIS
perl cdHitClstr2table.pl -s sampleList -c cd-hit.clstr -o outputTable.tsv

 Options:
   --help [-h]              Print this help
   --sample_list [-s]       List of samples. A text file with one sample per line
   --input_clstr [-i]       CD-HIT-EST ".clstr" output file
   --output_table [-o]      QIIME-ready output table

=head1 OPTIONS

=over 4

=item B<--help [-h]>

Print this help

=item B<--sample_list [-s]>

List of samples. A text file with one sample per line

=item B<--input_clstr [-i]>

CD-HIT-EST ".clstr" output file

=item B<--output_table [-o]>

QIIME-ready output table

=back

=head1 DESCRIPTION

B<cdHitClstr2table.pl> will convert the output from CD-HIT-EST to a table ready to input on QIIME.
A sample list is also required as input to have a complete output table.

=cut
