#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Bio::SeqIO;
use Bio::SeqIO::embl;


my $filename = $ARGV[0];
my $stream = Bio::SeqIO->new(-file => $filename, -format => 'EMBL');

my $num_seq = 0;
my %hla;
my %partial_count;
my %black_count;
my $psuedo_gene = 0;
while ( (my $seq = $stream->next_seq()) ) {
    # do something with $seq
    my @source = $seq->get_SeqFeatures('source');
    my @cds = $seq->get_SeqFeatures('CDS');
    if (scalar @cds == 1) {
        my $cds_tag = $cds[0];
        my @gene = $cds_tag->get_tag_values('gene');
        ++$hla{$gene[0]} if scalar @gene == 1;
        ++$partial_count{$gene[0]} if $cds_tag->has_tag('partial');
        if (scalar @source == 1) {
            my $source_tag = $source[0];
            if ($source_tag->has_tag('ethnic')) {
                my @ethnic = $source_tag->get_tag_values('ethnic');
                if (scalar @ethnic == 1) {
                    ++$black_count{$gene[0]} if $ethnic[0] =~ m/Black/;
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
print "of which ".$hla{$key}." are from $key\n";
print "of which ".$partial_count{$key}." are partial sequences\n" if exists $partial_count{$key};
print "of which ".$black_count{$key}." are sequences with Ethnic tag \"Black\"\n" if exists $black_count{$key};
}
print "Psuedo genes $psuedo_gene\n";
