#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Bio::SeqIO;
use Bio::SeqIO::embl;
use List::Util qw(reduce);
use Text::CSV;
 
my %gene_list;
my $csv = Text::CSV->new ( { sep_char => "\t" } )
                or die "Cannot use CSV: ".Text::CSV->error_diag ();
 
open my $fh, "<", "official_ref.txt" or die "official_ref.txt: $!";
while ( my $row = $csv->getline( $fh ) ) {
    $gene_list{"id"}{$row->[0]} = $row->[1];
    $gene_list{"acc"}{$row->[0]} = $row->[2];
}
$csv->eof or $csv->error_diag();
close $fh;
 
$csv->eol ("\r\n");

my $filename = $ARGV[0];
my $stream = Bio::SeqIO->new(-file => $filename, -format => 'EMBL');

my $num_seq = 0;
my %hla = ();
my %acc_idx = ();

while ( (my $seq = $stream->next_seq()) ) {
    # do something with $seq
    my @cds = $seq->get_SeqFeatures('CDS');
    if (scalar @cds == 1) {
        my $cds_tag = $cds[0];
        my @gene = $cds_tag->get_tag_values('gene');
        my @allele = $cds_tag->get_tag_values('allele');
        if (scalar @gene == 1 && ($gene[0] =~ /HLA-[A-C]/ || $gene[0] =~ /HLA-D[QPR].*/ )) {
            my $gene_name = $gene[0];
            my $allele_name = $allele[0];
            $hla{$gene_name}{$allele_name} = $seq->seq();
            $acc_idx{$seq->accession_number()} = $gene_name;
#            my @features = $seq->get_SeqFeatures();
#            my @allele_split = split /\*/, $allele_name;
#            my @allele_breakdown = split /:/, $allele_split[1];
#            foreach my $feature (@features) {
#                if ($feature->primary_tag eq  'exon') {
#                    my @number = $feature->get_tag_values('number');
#                    if ($number[0] eq '2' || $number[0] eq '3' || $number[0] eq '4') {
#                         $hla{$allele_name}{$number[0]} = $feature->seq()->seq();
#                    }
#                }
#            }
       }
    } else {
#        ++$psuedo_gene;
#It's a psuedo gene do nothing
    }
    ++$num_seq;
}

my $lef = $gene_list{"id"};
foreach my $key (sort keys %{$lef}) {
    print $key."\n";
    my @a = keys %{$hla{$key}};
    print $hla{$key}{$a[0]}."\n";
}

print STDERR "$num_seq processed\n";
