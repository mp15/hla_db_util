#!/usr/bin/env perl

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

foreach my $key (keys %ref_in) {
	my $value = $ref_in{$key};
	my @fields = split(/:/,$value);
	my $two = join(':', @fields[0..1]);
	$ref_in{$two} = $value if !exists($ref_in{$two});
}

my $input_fh;
open($input_fh, "<-") or die "Cannot open stdin";
while (<$input_fh>) {
	chomp;
	if (exists($ref_in{$_})) {
		print $ref_in{$_}."\n";
	} else {
		print $_."\n";
		print STDERR "Warning unknown allele: $_\n"
	}
}
close $input_fh;
