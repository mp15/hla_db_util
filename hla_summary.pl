#!/usr/bin/env perl

#  hla_summary.pl
#  hla_db_util - set of utilities for working with IMGT/HLA database data
#
# Created by Martin Pollard.
# Copyright © 2015, 2016, 2017 Genome Research Limited.
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
print STDERR "There are $num_seq sequences in this file\n";

print "\t";
foreach my $key (sort keys %hla) {
    print "\t${key}";
}
print "\n";
foreach my $count (2..4) {
    print "partial\t";
    foreach my $key (sort keys %hla) {
        if (exists $partial_count{$key}) {
            print "\t".$partial_count{$key};
        } else {
            print "\t0";
        }
    }
    print "\n";
    print "total\t";
    foreach my $key (sort keys %hla) {
        if (exists $hla{$key}) {
            print "\t".$hla{$key};
        } else {
            print "\t0";
        }
    }
    print "\n";
}

print STDERR "Psuedo genes $psuedo_gene\n";
