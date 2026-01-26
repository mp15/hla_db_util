#!/usr/bin/env perl

#  hla_exons_fasta.pl
#  hla_db_util - set of utilities for working with IMGT/HLA database data
#
# Created by Martin Pollard.
# Copyright © 2018 Genome Research Limited.
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
        my @allele = $cds_tag->get_tag_values('allele');
        if (scalar @gene == 1 && ($gene[0] =~ /HLA-[A-C]/ || $gene[0] =~ /HLA-D[QPR].*/ )) {
            my $gene_name = $gene[0];
            my $allele_name = $allele[0];
            my @features = $seq->get_SeqFeatures();
            my @allele_split = split /\*/, $allele_name;
            my @allele_breakdown = split /:/, $allele_split[1];
            foreach my $feature (@features) {
                if ($feature->primary_tag eq  'exon') {
                    my @number = $feature->get_tag_values('number');
                    if ($number[0] eq '2' || $number[0] eq '3' || $number[0] eq '4') {
                         $hla{$allele_name}{$number[0]} = $feature->seq()->seq();
                    }
                }
            }
       }
    } else {
        ++$psuedo_gene;
    }
    ++$num_seq;
}
print STDERR "There are $num_seq sequences in this file\n";
foreach my $key (sort keys %hla) {
next if ($key !~ /^HLA-[A-C].*/);
print ">${key}\n";
print $hla{$key}{2};
for (my $i=0; $i<100;++$i) {print 'N';}
print $hla{$key}{3};
for (my $i=0; $i<100;++$i) {print 'N';}
print $hla{$key}{4} if exists $hla{$key}{4};
print "\n";
}
foreach my $key (sort keys %hla) {
next if ($key !~ /^HLA-D(RB[1345]|[PQ][AB]1).*/);
next if ($key =~ /^HLA-DRB4\*03:01N/); # No exon 2 data
print ">${key}\n";
print $hla{$key}{2} if exists $hla{$key}{2};
for (my $i=0; $i<100;++$i) {print 'N';}
print $hla{$key}{3} if exists $hla{$key}{3};
print"\n";
}
print STDERR "Psuedo genes $psuedo_gene\n";
