#!/usr/bin/env perl

#  g_group.pl
#  hla_db_util - set of utilities for working with IMGT/HLA database data
#
# Created by Martin Pollard.
# Copyright © 2016, 2017 Genome Research Limited.
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

# Usage g_group.pl <ref_g_groups>
# reads input from STDIN and writes converted output to STDOUT

my $ref = $ARGV[0];

my $ref_fh;
my %ref_in;
open($ref_fh, "<", $ref) or die "Could not open reference $ref";

while (<$ref_fh>) {
	chomp;
	next if $_ =~ /^#/; #skip header lines and initial comments
	my @a = split /;/;
	if (scalar @a == 2) { # 2 fields means there is no G group
		my @fields = split(/:/,$a[1]);
		if (scalar @fields == 4) { # if this is a four field allele need to create a g group at 3 fields
			my $three = join(':', @fields[0..2]);
			$ref_in{"HLA-$a[0]$a[1]"} = "HLA-$a[0]${three}G";
		} else {
			$ref_in{"HLA-$a[0]$a[1]"} = "HLA-$a[0]$a[1]";
		}
	} else {
		foreach my $b (split(/\//, $a[1])) {
			my @fields = split(/:/,$b);
			if (scalar @fields == 4) {
				my $three = join(':', @fields[0..2]);
				$ref_in{"HLA-$a[0]$three"} = "HLA-$a[0]$a[2]";
			}

			$ref_in{"HLA-$a[0]$b"} = "HLA-$a[0]$a[2]";
		}
	}
}
close $ref_fh;

foreach my $key (sort keys %ref_in) {
	my $value = $ref_in{$key};
	my @fields = split(/:/,$value);
	my $two = join(':', @fields[0..1]);
	$ref_in{$two} = $value if !exists($ref_in{$two});
}

sub translate($$)
{
    my ($ref_data, $query) = @_;
    if (exists($ref_data->{$query})) {
        return($ref_data->{$query});
    } else {
        print STDERR "Warning unknown allele: ${query}\n";
        return($query);
    }
}

open(my $input_fh, "<-") or die "Cannot open stdin";
while (<$input_fh>) {
	chomp;
	# Pass through missing allele indicators
	if ( $_ eq "****" ) {
		print "****\n";
		next;
	}

    print translate(\%ref_in,$_)."\n";
}
close $input_fh;
