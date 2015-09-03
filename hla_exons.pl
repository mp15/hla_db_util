#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
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
        if (scalar @gene == 1) {
            my $gene_name = $gene[0];
            my @features = $seq->get_SeqFeatures();
            my $accum = "";
            foreach my $feature (@features) {
                if ($feature->primary_tag eq 'UTR') {
                    $accum = "${accum}U";
                } elsif ($feature->primary_tag eq  'exon') {
                    my @number = $feature->get_tag_values('number');
                    $accum = "${accum}".$number[0];
                } elsif ($feature->primary_tag eq 'intron') {
                    $accum = "${accum}I"; 
                }
            }
            $hla{$gene_name}{$accum} += 1;

#            my @exons = map {$_->get_tag_values('number')} $seq->get_SeqFeatures('exon');
#            my @exon_bits = map { 1<<($_-1)} @exons;
#            my $exon_bitmask = reduce {$a |= $b} @exon_bits;
#            $hla{$gene_name}{$exon_bitmask} += 1;
#            my @utr_count = $seq->get_SeqFeatures('UTR');
#            $utr{$gene_name}{scalar @utr_count} += 1;
        }
    } else {
        ++$psuedo_gene;
    }
    ++$num_seq;
}
print "There are $num_seq sequences in this file\n";
foreach my $key (sort keys %hla) {
print "[$key]\n";
foreach my $allele (sort { length $a <=> length $b } keys $hla{$key}) {
    print "$allele = ".$hla{$key}{$allele}."\n";
}
#foreach my $vec (sort keys $hla{$key}) {
#my $exons = sprintf("%11b", $vec);
#$exons = reverse $exons;
#print "$exons = ".$hla{$key}{$vec}."\n";
#}
#print "UTR\n";
#foreach my $count (sort keys $utr{$key}) {
#print "$count = ".$utr{$key}{$count}."\n";
#}
}
print "Psuedo genes $psuedo_gene\n";
