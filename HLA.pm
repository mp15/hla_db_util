#!/usr/bin/env perl

#  HLA.pm
#  hla_db_util - set of utilities for working with IMGT/HLA database data
#
# Created by Martin Pollard.
# Copyright © 2016, 2017, 2026 Genome Research Limited.
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


package HLA;
use strict;
use warnings;
use B qw(svref_2object);


use Exporter qw(import);

our @EXPORT_OK = qw(print_full_matches print_three_field_matches print_two_field_matches print_one_field_matches parse_soap parse_flat parse_file parse_multiresfile print_full_matches_multi print_one_field_matches_multi print_two_field_matches_multi print_three_field_matches_multi);

sub three_field_compare($$) {
    my ($allele_a, $allele_b) = @_;

    my @allele_a_gene_allele = split /\*/, $allele_a;
    my @allele_b_gene_allele = split /\*/, $allele_b;
    return 0 if (scalar @allele_a_gene_allele != 2 || scalar @allele_b_gene_allele != 2);
    return 0 if ($allele_a_gene_allele[0] ne $allele_b_gene_allele[0]);

    my @a_fields = split /:/, $allele_a_gene_allele[1];
    my @b_fields = split /:/, $allele_b_gene_allele[1];
    return 1 if (scalar @a_fields == 1) && (scalar @b_fields == 1) && ($a_fields[0] eq $b_fields[0]);
    return 1 if (scalar @a_fields == 2) && (scalar @b_fields == 2) && ($a_fields[0] eq $b_fields[0]) && ($a_fields[1] eq $b_fields[1]);
    return 0 if ((scalar @a_fields < 3)
        || (scalar @b_fields < 3)
        || ($a_fields[0] ne $b_fields[0])
        || ($a_fields[1] ne $b_fields[1])
        || ($a_fields[2] ne $b_fields[2]));

    return 1;
}

sub two_field_compare($$) {
	my ($allele_a, $allele_b) = @_;

	my @allele_a_gene_allele = split /\*/, $allele_a;
	my @allele_b_gene_allele = split /\*/, $allele_b;
	return 0 if (scalar @allele_a_gene_allele != 2 || scalar @allele_b_gene_allele != 2);
	return 0 if ($allele_a_gene_allele[0] ne $allele_b_gene_allele[0]);

	my @a_fields = split /:/, $allele_a_gene_allele[1];
	my @b_fields = split /:/, $allele_b_gene_allele[1];
    return 1 if (scalar @a_fields == 1) && (scalar @b_fields == 1) && ($a_fields[0] eq $b_fields[0]);
	return 0 if ((scalar @a_fields < 2)
        || (scalar @b_fields < 2)
        || ($a_fields[0] ne $b_fields[0])
        || ($a_fields[1] ne $b_fields[1]));

	return 1;
}

sub one_field_compare($$) {
	my ($allele_a, $allele_b) = @_;

	my @allele_a_gene_allele = split /\*/, $allele_a;
	my @allele_b_gene_allele = split /\*/, $allele_b;
	return 0 if (scalar @allele_a_gene_allele != 2 || scalar @allele_b_gene_allele != 2);
	return 0 if ($allele_a_gene_allele[0] ne $allele_b_gene_allele[0]);

	my @a_fields = split /:/, $allele_a_gene_allele[1];
	my @b_fields = split /:/, $allele_b_gene_allele[1];
	return 0 if ($a_fields[0] ne $b_fields[0]);

	return 1;
}

sub full_field_compare($$) {
	my ($allele_a, $allele_b) = @_;
	my (@allele_a_gene_allele) = split /\//, $allele_a;
	my (@allele_b_gene_allele) = split /\//, $allele_b;
	my %count = ();

	foreach my $element (@allele_a_gene_allele, @allele_b_gene_allele) { $count{$element}++ }
	foreach my $element (keys %count) {
		if ($count{$element} > 1) {
			return 1;
		}
	}

	return 0;
}

sub print_full_matches($\%\%) {
	my $compare_func = \&full_field_compare;
	unshift(@_, $compare_func);
	print_n_field_matches(@_);
}

sub print_three_field_matches($\%\%) {
	my $compare_func = \&three_field_compare;
	unshift(@_, $compare_func);
	print_n_field_matches(@_);
}

sub print_two_field_matches($\%\%) {
	my $compare_func = \&two_field_compare;
	unshift(@_, $compare_func);
	print_n_field_matches(@_);
}

sub print_one_field_matches($\%\%) {
	my $compare_func = \&one_field_compare;
	unshift(@_, $compare_func);
	print_n_field_matches(@_);
}

sub print_n_field_matches($$\%\%) {
	my ($compare_func,$output_n_name,$sanger_data_h, $hiseqx_data_h) = @_;
	open (my $output_n_fh, ">", $output_n_name) or die "Cannot open output file for ".svref_2object($compare_func)->GV->NAME." matches $output_n_name";
	my %allele_count;
	my $allele;
	my $sample;
	my %sanger_data = %$sanger_data_h;
	my %hiseqx_data = %$hiseqx_data_h;

	foreach $allele (keys %hiseqx_data) {
		$allele_count{$allele}{matches} = 0;
		$allele_count{$allele}{total} = 0;
                $allele_count{$allele}{missing_sanger} = 0;
                $allele_count{$allele}{missing_other} = 0;
		next if (!exists $sanger_data{$allele} || !exists $hiseqx_data{$allele});
		foreach $sample (keys %{$hiseqx_data{$allele}}) {
			if (exists $sanger_data{$allele}{$sample} && exists $hiseqx_data{$allele}{$sample}) {
				if (scalar @{$sanger_data{$allele}{$sample}} == 2 && scalar @{$hiseqx_data{$allele}{$sample}} == 2) {
					if ((&$compare_func($sanger_data{$allele}{$sample}[0],$hiseqx_data{$allele}{$sample}[0]) &&
						&$compare_func($sanger_data{$allele}{$sample}[1], $hiseqx_data{$allele}{$sample}[1])) ||
						(&$compare_func($sanger_data{$allele}{$sample}[0], $hiseqx_data{$allele}{$sample}[1]) &&
						&$compare_func($sanger_data{$allele}{$sample}[1], $hiseqx_data{$allele}{$sample}[0]))) {
						$allele_count{$allele}{matches} += 2;
					} elsif ((&$compare_func($sanger_data{$allele}{$sample}[0], $hiseqx_data{$allele}{$sample}[0]) ||
						&$compare_func($sanger_data{$allele}{$sample}[1], $hiseqx_data{$allele}{$sample}[1])) ||
						(&$compare_func($sanger_data{$allele}{$sample}[0], $hiseqx_data{$allele}{$sample}[1]) ||
						&$compare_func($sanger_data{$allele}{$sample}[1], $hiseqx_data{$allele}{$sample}[0]))) {
						$allele_count{$allele}{matches} += 1;
					} else { print "mismatch2of2\t".svref_2object($compare_func)->GV->NAME."\t".$sample."\t".$sanger_data{$allele}{$sample}[0]."\t".$sanger_data{$allele}{$sample}[1],"\t".$hiseqx_data{$allele}{$sample}[0]."\t".$hiseqx_data{$allele}{$sample}[1]."\n";}
				} elsif (scalar @{$sanger_data{$allele}{$sample}} == 1 && scalar @{$hiseqx_data{$allele}{$sample}} == 2) {
					if (&$compare_func($sanger_data{$allele}{$sample}[0],$hiseqx_data{$allele}{$sample}[0]) &&
						&$compare_func($sanger_data{$allele}{$sample}[0], $hiseqx_data{$allele}{$sample}[1])) {
						$allele_count{$allele}{matches} += 2;
					} elsif (&$compare_func($sanger_data{$allele}{$sample}[0],$hiseqx_data{$allele}{$sample}[0]) ||
						&$compare_func($sanger_data{$allele}{$sample}[0], $hiseqx_data{$allele}{$sample}[1])) {
						$allele_count{$allele}{matches} += 1;
					} else {}
				} elsif (scalar @{$sanger_data{$allele}{$sample}} == 2 && scalar @{$hiseqx_data{$allele}{$sample}} == 1) {
					if (&$compare_func($sanger_data{$allele}{$sample}[0],$hiseqx_data{$allele}{$sample}[0]) &&
						&$compare_func($sanger_data{$allele}{$sample}[1], $hiseqx_data{$allele}{$sample}[0])) {
						$allele_count{$allele}{matches} += 2;
					} elsif (&$compare_func($sanger_data{$allele}{$sample}[0],$hiseqx_data{$allele}{$sample}[0]) ||
						&$compare_func($sanger_data{$allele}{$sample}[1], $hiseqx_data{$allele}{$sample}[0])) {
						$allele_count{$allele}{matches} += 1;
					} else {}
				} elsif (scalar @{$sanger_data{$allele}{$sample}} == 1 && scalar @{$hiseqx_data{$allele}{$sample}} == 1) {
					if (&$compare_func($sanger_data{$allele}{$sample}[0],$hiseqx_data{$allele}{$sample}[0])) {
						$allele_count{$allele}{matches} += 2;
					} else {}
				}

				$allele_count{$allele}{total} += 2;
			} elsif (exists $sanger_data{$allele}{$sample} ) {
                                $allele_count{$allele}{missing_sanger} += 2;
                        } elsif (exists $hiseqx_data{$allele}{$sample}) {
                                $allele_count{$allele}{missing_other} += 2;
                        }
		}
	}

	print $output_n_fh "Gene\tMatches\tTotal\tProportion\tMissing SBT\tMissing Compared\n";
	foreach $allele (sort keys %hiseqx_data) {
		next if (!exists $sanger_data{$allele}) or ($allele_count{$allele}{total} == 0);
		printf($output_n_fh "${allele}\t%d\t%d\t%.3f\t%d\t%d\n",$allele_count{$allele}{matches},$allele_count{$allele}{total},$allele_count{$allele}{matches}/$allele_count{$allele}{total},$allele_count{$allele}{missing_sanger},$allele_count{$allele}{missing_other});
	}
	close($output_n_fh);
}

sub parse_soap($) {
	my $soap = shift;
	open (my $soap_fh, "<", $soap) or die "Cannot open SOAP data $$soap";
	my %soap_data;
	while (my $input = <$soap_fh>) {
		chomp $input;
		my @fields = split /\t/, $input;
		next if $fields[2] eq '---'; #handle the null case
		$fields[1] = "HLA-$fields[1]";
		$fields[2] = "HLA-$fields[2]";
		$fields[3] = "HLA-$fields[3]";
		if (exists $soap_data{$fields[1]} && exists $soap_data{$fields[1]}{$fields[0]}) {
			push @{$soap_data{$fields[1]}{$fields[0]}}, $fields[2];
		} else {
			$soap_data{$fields[1]}{$fields[0]} = [$fields[2]];
		}
	}
	close($soap_fh);
	return %soap_data;
}

sub parse_flat($) {
	my $flat = shift;
	open (my $flat_fh, "<", $flat) or die "Cannot open flat data $$flat";
	my %flat_data;
	while (my $input = <$flat_fh>) {
		chomp $input;
		my @fields = split /\t/, $input;
		next if $fields[2] eq '0000'; #handle the null case
		$fields[1] = "$fields[1]";
		$fields[2] = "$fields[2]";
		if (exists $flat_data{$fields[1]} && exists $flat_data{$fields[1]}{$fields[0]}) {
			push @{$flat_data{$fields[1]}{$fields[0]}}, $fields[2];
		} else {
			$flat_data{$fields[1]}{$fields[0]} = [$fields[2]];
		}
	}
	close($flat_fh);
	return %flat_data;
}

sub parse_file($) {
	my $sanger=shift;
	open(my $sangerfh, "<", $sanger) or die "can't open $sanger";
	# header
	my $raw_header = <$sangerfh>;
	chomp $raw_header;
	my @sanger_header = split /\t/, $raw_header;

	my %sanger_data;

	while (my $input = <$sangerfh>) {
		chomp $input;
		my @fields = split /\t/, $input;
		for (my $i = 1; $i < scalar @fields; ++$i) {
			next if $fields[$i] =~ '-$';
			if (exists $sanger_data{$sanger_header[$i]} && exists $sanger_data{$sanger_header[$i]}{$fields[0]}) {
				push @{$sanger_data{$sanger_header[$i]}{$fields[0]}}, $fields[$i];
			} else {
				$sanger_data{$sanger_header[$i]}{$fields[0]} = [$fields[$i]];
			}
		}
	}

	close $sangerfh;
	return %sanger_data;
}

sub parse_multiresfile($) {
	my $sanger=shift;
	open(my $sangerfh, "<", $sanger) or die "can't open $sanger";
	# first three lines are header then blank line so we'll just read until blank line
	while(($_ = <$sangerfh>) =~ /\S/) {}

	my %sanger_data;
	my $current_sample;
	while (my $input = <$sangerfh>) {
		chomp $input;
		my @fields = split /\t/, $input;
		next if scalar @fields == 1;
		chomp $fields[0];
		$current_sample = $fields[0] if $fields[0] !~ /^\s*$/;

		next if $fields[1] =~ 'NA$';
		my @allele = split /\*/, $fields[1];
		my $gene = "HLA-$allele[0]";
		my @pair = ("HLA-$fields[1]","HLA-$fields[2]");
		if (exists $sanger_data{$gene} && exists $sanger_data{$gene}{$current_sample}) {
			push @{$sanger_data{$gene}{$current_sample}}, \@pair;
		} else {
			$sanger_data{$gene}{$current_sample} = [\@pair];
		}
	}

	close $sangerfh;
	return %sanger_data;
}

sub print_full_matches_multi($$\%\%) {
	my $compare_func = \&full_field_compare;
	unshift(@_, $compare_func);
	print_n_field_matches_multi(@_);
}

sub print_three_field_matches_multi($$\%\%) {
    my $compare_func = \&three_field_compare;
    unshift(@_, $compare_func);
    print_n_field_matches_multi(@_);
}

sub print_two_field_matches_multi($$\%\%) {
	my $compare_func = \&two_field_compare;
	unshift(@_, $compare_func);
	print_n_field_matches_multi(@_);
}

sub print_one_field_matches_multi($$\%\%) {
	my $compare_func = \&one_field_compare;
	unshift(@_, $compare_func);
	print_n_field_matches_multi(@_);
}

sub print_n_field_matches_multi($$\%\%) {
	my ($compare_func,$output_two_name,$error_output,$sanger_data_h, $hiseqx_data_h) = @_;
	open (my $output_full_fh, ">", $output_two_name) or die "Cannot open output file for full matches $output_two_name";
	my %allele_count;
	my $allele;
	my $sample;
	my %sanger_data = %$sanger_data_h;
	my %hiseqx_data = %$hiseqx_data_h;

	foreach $allele (keys %hiseqx_data) {
		$allele_count{$allele}{matches} = 0;
		$allele_count{$allele}{total} = 0;
                $allele_count{$allele}{missing_sbt} = 0;
                $allele_count{$allele}{missing_other} = 0;
                $allele_count{$allele}{missing_both} = 0;
		next if (!exists $sanger_data{$allele} || !exists $hiseqx_data{$allele});
		foreach $sample (keys %{$hiseqx_data{'HLA-A'}}) {
			if (exists $sanger_data{$allele}{$sample} && exists $hiseqx_data{$allele}{$sample}) {
				my @fixed_hiseqx;
				if (scalar @{$hiseqx_data{$allele}{$sample}} == 1) {
					@fixed_hiseqx = ($hiseqx_data{$allele}{$sample}[0], $hiseqx_data{$allele}{$sample}[0]);
				} else {
					@fixed_hiseqx = ($hiseqx_data{$allele}{$sample}[0], $hiseqx_data{$allele}{$sample}[1]);
				}
				my $match_found = 0;
				foreach my $pairref (@{$sanger_data{$allele}{$sample}}) {
					if ((&$compare_func(@{$pairref}[0],$fixed_hiseqx[0]) &&
						&$compare_func(@{$pairref}[1], $fixed_hiseqx[1])) ||
						(&$compare_func(@{$pairref}[0], $fixed_hiseqx[1]) &&
						&$compare_func(@{$pairref}[1], $fixed_hiseqx[0]))) {
						$allele_count{$allele}{matches} += 2;
						$match_found = 1;
						last;
					}
				}
				print $error_output "match not found\t".svref_2object($compare_func)->GV->NAME."\t$sample\t$allele\t$fixed_hiseqx[0]\t$fixed_hiseqx[1]\n" if $match_found == 0;
				$allele_count{$allele}{total} += 2;
			} elsif (exists $sanger_data{$allele}{$sample} ) {
                                print $error_output "missing_other\t".svref_2object($compare_func)->GV->NAME."\t$sample\t$allele\t-\t-\n";
				$allele_count{$allele}{missing_other} += 2;
                        } elsif (exists $hiseqx_data{$allele}{$sample}) {
                                print $error_output "missing_sbt\t".svref_2object($compare_func)->GV->NAME."\t$sample\t$allele\t-\t-\n";
                                $allele_count{$allele}{missing_sbt} += 2;
                        } else {
				print $error_output "missing_both\t".svref_2object($compare_func)->GV->NAME."\t$sample\t$allele\t-\t-\n";
				$allele_count{$allele}{missing_both} += 2;
			}
		}
	}

	print $output_full_fh "Gene\tMatches\tTotal\tProportion\tMissing SBT\tMissing Other\tMissing Both\n";
	foreach $allele (sort keys %hiseqx_data) {
		next if (!exists $sanger_data{$allele});
		printf($output_full_fh "${allele}\t".$allele_count{$allele}{matches}."\t".$allele_count{$allele}{total}."\t"."%.3f\t".$allele_count{$allele}{missing_sbt}."\t".$allele_count{$allele}{missing_other}."\t".$allele_count{$allele}{missing_both}."\n",$allele_count{$allele}{matches}/$allele_count{$allele}{total});
	}
	close($output_full_fh);
}

1;
