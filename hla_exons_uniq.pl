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
my %names = ();
my %hla_n = ();
my %names_n = ();
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
            if ( $gene[0] =~ /HLA-[A-C]/ ) {
                my $exon2;
                my $exon3;
                my $exon2_n;
                my $exon3_n;
                foreach my $feature (@features) {
                    if ($feature->primary_tag eq  'exon') {
                        my @number = $feature->get_tag_values('number');
                        if ($number[0] eq '2') {
                            #hold it!
                            $exon2 = $feature->seq()->translate->seq();
                            $exon2_n = $feature->seq()->seq();
                        } elsif ($number[0] eq '3') {
                            $exon3 = $feature->seq()->translate->seq();
                            $exon3_n = $feature->seq()->seq();
                        }
                    }
                }
                if ($exon2 && $exon3) {
                    $hla{$gene_name}{$exon2."Z".$exon3}++;
                    $names{$gene_name}{$exon2."Z".$exon3}{$allele_breakdown[0].":".$allele_breakdown[1]}++;
                    $hla_n{$gene_name}{$exon2_n."Z".$exon3_n}++;
                    $names_n{$gene_name}{$exon2_n."Z".$exon3_n}{$allele_breakdown[0].":".$allele_breakdown[1]}++;
                }
            } else { # $gene[0] =~ /HLA-D[QPR].*/
                my $exon2;
                my $exon2_n;
                foreach my $feature (@features) {
                    if ($feature->primary_tag eq  'exon') {
                        my @number = $feature->get_tag_values('number');
                        if ($number[0] eq '2') {
                            #hold it!
                            $exon2 = $feature->seq()->translate->seq();
                            $exon2_n = $feature->seq()->seq();
                        }
                    }
                }
                if ($exon2) {
                    $hla{$gene_name}{$exon2}++;
                    $names{$gene_name}{$exon2}{$allele_breakdown[0].":".$allele_breakdown[1]}++;
                    $hla_n{$gene_name}{$exon2_n}++;
                    $names_n{$gene_name}{$exon2_n}{$allele_breakdown[0].":".$allele_breakdown[1]}++;
                }
            }
       }
    } else {
        ++$psuedo_gene;
    }
    ++$num_seq;
}
print STDERR "There are $num_seq sequences in this file\n";
my %transpose;
foreach my $key (sort keys %hla) {
    next if ($key !~ /^HLA-[A-C].*/ && $key !~ /HLA-D[QPR].*/);
    print "$key\n";
    my $unique = 0;
    my $nonunique = 0;
    foreach my $seq (sort keys $hla{$key}) {
        if ($hla{$key}{$seq} == 1) {
            ++$unique;
            next;
        }
        ++$nonunique;
        #print $hla{$key}{$seq}."\n";
        #foreach my $allele (sort keys $names{$key}{$seq}) {
        #    print $allele." ".$names{$key}{$seq}{$allele}."\n";
        #}
    }
    print "nonunique: $nonunique\n";
    print "unique: $unique\n";
    $transpose{$key}{"unique"} = $unique;
    $transpose{$key}{"nonunique"} = $nonunique;
}
my %transpose_n;
foreach my $key (sort keys %hla_n) {
    next if ($key !~ /^HLA-[A-C].*/ && $key !~ /HLA-D[QPR].*/);
    print "$key\n";
    my $unique = 0;
    my $nonunique = 0;
    foreach my $seq (sort keys $hla_n{$key}) {
        if ($hla_n{$key}{$seq} == 1) {
            ++$unique;
            next;
        }
        ++$nonunique;
        #print $hla{$key}{$seq}."\n";
        #foreach my $allele (sort keys $names{$key}{$seq}) {
        #    print $allele." ".$names{$key}{$seq}{$allele}."\n";
        #}
    }
    print "nonunique: $nonunique\n";
    print "unique: $unique\n";
    $transpose_n{$key}{"unique"} = $unique;
    $transpose_n{$key}{"nonunique"} = $nonunique;
}

print "protein\n";
print "Allele";
foreach my $allele (sort keys %transpose) {
    print "\t$allele";
}
print "\n";
print "Nonunique";
foreach my $allele (sort keys %transpose) {
    print "\t".$transpose{$allele}{"nonunique"};
}
print "\n";
print "Unique";
foreach my $allele (sort keys %transpose) {
    print "\t".$transpose{$allele}{"unique"};
}
print "\n";

print "nuc\n";
print "Allele";
foreach my $allele (sort keys %transpose_n) {
    print "\t$allele";
}
print "\n";
print "Nonunique";
foreach my $allele (sort keys %transpose_n) {
    print "\t".$transpose_n{$allele}{"nonunique"};
}
print "\n";
print "Unique";
foreach my $allele (sort keys %transpose_n) {
    print "\t".$transpose_n{$allele}{"unique"};
}
print "\n";

    #foreach my $key (sort keys %hla) {
    #next if ($key !~ /^HLA-D(RB[1345]|[PQ][AB]1).*/);
    #next if ($key =~ /^HLA-DRB4*03:01N/; # No exon 2 data
    #}
print STDERR "Psuedo genes $psuedo_gene\n";
