#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper qw(Dumper);
use Bio::SeqIO;
use Bio::SeqIO::embl;
use List::Util qw(reduce);

 
my $filename = $ARGV[0];
my $stream = Bio::SeqIO->new(-file => $filename, -format => 'EMBL');

my $num_seq = 0;
my %hla = ();
my %utr = ();
my %partial_count;
my $psuedo_gene = 0;
while ( (my $seq = $stream->next_seq()) ) {
    # do something with $seq
    my @cds = $seq->get_SeqFeatures('CDS');
    if (scalar @cds == 1) {
        my $cds_tag = $cds[0];
        my @gene = $cds_tag->get_tag_values('gene');
        my @allele = $cds_tag->get_tag_values('allele');
        if (scalar @gene == 1 && ($gene[0] =~ /HLA-[A-C]/ || $gene[0] =~ /HLA-D[QPR][AB][1-5].*/ )) {
            my $gene_name = $gene[0];
            my $allele_name = $allele[0];
            my @features = $seq->get_SeqFeatures();
            my @allele_split = split /\*/, $allele_name;
            my @allele_breakdown = split /:/, $allele_split[1];
            ++$hla{$gene_name}{scalar @allele_breakdown};
       }
    } else {
        ++$psuedo_gene;
    }
    ++$num_seq;
}
print "There are $num_seq sequences in this file\n";
print "Digits/Alleles\t";
foreach my $key (sort keys %hla) {
    print "\t${key}";
}
print "\n";
foreach my $count (2..4) {
    print "$count\t";
    foreach my $key (sort keys %hla) {
        if (exists $hla{$key}{$count}) {
            print "\t".$hla{$key}{$count};
        } else {
            print "\t0";
        }
    }
    print "\n";
}
print "Psuedo genes $psuedo_gene\n";
