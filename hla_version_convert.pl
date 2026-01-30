#!/usr/bin/env perl

#  hla_version_convert.pl
#  hla_db_util - set of utilities for working with IMGT/HLA database data
#
# Created by Martin Pollard.
# Copyright © 2016, 2018 Genome Research Limited.
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

# Convert HLA allele names from one version to another using a reference allele history file

use strict;
use warnings;

# Usage g_group.pl <ref_Allele_history> <from version> <to version>
# reads input from STDIN and writes converted output to STDOUT

my $ref = $ARGV[0];
my $from = $ARGV[1];
my $to = $ARGV[2];
my $hard = $ARGV[3];
my $verbose = $ARGV[4];

my $ref_fh;
my %ref_in;
open($ref_fh, "<", $ref) or die "Could not open allele_history $ref";

# Read header line
my @header;
while (<$ref_fh>) {
	chomp;

	# Find final header line
	if (/^HLA_ID/) {
		@header = split(/,/,$_);
		last;
	}
}
if (scalar(@header) == 0) {
	die "Could not find header line in allele history file $ref\n";
}

my $i = 0;
my %version_id;
foreach my $version (@header) {
	$version_id{$version} = $i++;
}
if (!exists($version_id{$from})) {
	die "Unknown from version $from. It should be input without the '.'s.\n";
}
if (!exists($version_id{$to})) {
	die "Unknown to version $to. It should be input without the '.'s.\n";
}
my $from_id = $version_id{$from};
my $to_id = $version_id{$to};
my %translate_to_hlaid;
my %hlaid_to_hlaid;
my %hlaid_to_translation;
while (<$ref_fh>) {
	chomp;
	my @a = split /,/;
	$translate_to_hlaid{"HLA-$a[$from_id]"} = "$a[0]";
	$hlaid_to_translation{"$a[0]"} = "HLA-$a[$to_id]";
}
close $ref_fh;
#manual hlaid to hlaid translations
$hlaid_to_hlaid{"HLA00407"} = "HLA02120"; # Sequence shown to be identical to C*02:10:01:01 (May 2006)
$hlaid_to_hlaid{"HLA00481"} = "HLA04311"; # Sequence has now been shown to be in error and is identical to C*17:01:01:02 (August 2016)
$hlaid_to_hlaid{"HLA00332"} = "HLA14088"; #Sequence shown to contain errors and be identical to B*47:01:01:03 (April 2018)
$hlaid_to_hlaid{"HLA00624"} = "HLA29352"; # Sequence extended and renamed DQB1*02:180:01 (November 2020)

# Process input alleles
my $input_fh;
open($input_fh, "<-") or die "Cannot open stdin";
while (<$input_fh>) {
    chomp;
	# Pass through missing allele indicators
	if ( $_ eq "****" ) {
		print "****\n";
		next;
	}
	
	# Translate allele to HLA ID and then to target version, checking for overrides in between
    if (exists($translate_to_hlaid{$_})) {
		my $translated_hlaid = $translate_to_hlaid{$_};
		# Check for manual translation
		if (exists($hlaid_to_hlaid{$translate_to_hlaid{$_}})) {
			$translated_hlaid = $hlaid_to_hlaid{$translate_to_hlaid{$_}};
		}
		if (!exists($hlaid_to_translation{$translated_hlaid})) {
			die "Warning no translation for allele $_ ($translated_hlaid)\n" if $hard;
			print STDERR "N ($_ $translated_hlaid)\n" if $verbose;
		} else {
			print "$hlaid_to_translation{$translated_hlaid}\n";
			print STDERR "T ($_)\n" if $verbose;
		}
    } else {
		die "Warning unknown allele: $_\n" if $hard;
		print STDERR "U$_\n" if $verbose;
    }
}
close $input_fh;