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
        if (scalar @gene == 1 && ($gene[0] =~ /HLA-[A-C]/ || $gene[0] =~ /HLA-D[QPR].*/ )) {
            my $gene_name = $gene[0];
            my $allele_name = $allele[0];
            my @features = $seq->get_SeqFeatures();
            my @allele_split = split /\*/, $allele_name;
            my @allele_breakdown = split /:/, $allele_split[1];
++$hla{$gene_name}{scalar @allele_breakdown};
            foreach my $feature (@features) {
                if ($feature->primary_tag eq  'exon') {
                    my @number = $feature->get_tag_values('number');
                    if ($number[0] eq '2') {
#                         print ">".@allele_split[1]."\n".$feature->seq()->seq()."\n";
                    }
                }
            }
       }
    } else {
        ++$psuedo_gene;
    }
    ++$num_seq;
}
print "There are $num_seq sequences in this file\n";
foreach my $key (sort keys %hla) {
print "${key}\n";
foreach my $count (sort keys $hla{$key}) {
    print "$count = ".$hla{$key}{$count}."\n";
}
}
print "Psuedo genes $psuedo_gene\n";
