#!/usr/bin/env perl

#  hla_exons.pl
#  hla_db_util - set of utilities for working with IMGT/HLA database data
#
# Created by Martin Pollard.
# Copyright © 2015, 2018 Genome Research Limited.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


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
foreach my $allele (sort { length $a <=> length $b } keys %{$hla{$key}}) {
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
