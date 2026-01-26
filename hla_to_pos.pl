#!/usr/bin/env perl

#  hla_to_pos.pl
#  hla_db_util - set of utilities for working with IMGT/HLA database data
#
# Created by Martin Pollard.
# Copyright © 2016 Genome Research Limited.
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
