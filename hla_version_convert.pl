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

my $ref_fh;
my %ref_in;
open($ref_fh, "<", $ref) or die "Could not open allele_history $ref";


my @header;
while (<$ref_fh>) {
	chomp;

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
my %translate;
while (<$ref_fh>) {
	my @a = split /,/;
	$translate{"HLA-$a[$from_id]"} = "HLA-$a[$to_id]";
}

my $input_fh;
open($input_fh, "<-") or die "Cannot open stdin";
while (<$input_fh>) {
    chomp;
    if (exists($translate{$_})) {
        print "$translate{$_}\n";
	print STDERR "T ($_)\n";
    } else {
	die "Warning unknown allele: $_\n" if $hard;
	print STDERR "U$_\n";
    }
}
close $input_fh;

