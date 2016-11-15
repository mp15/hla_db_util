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
    next if $_ =~ /^#/;
	my @a = split /;/;
    if (scalar @a == 2) {
        $ref_in{"HLA-$a[0]$a[1]"} = "HLA-$a[0]$a[1]";
    } else {
        foreach my $b (split(/\//, $a[1])) {
            $ref_in{"HLA-$a[0]$b"} = "HLA-$a[0]$a[2]";
        }
    }
}
close $ref_fh;

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
