#!/usr/bin/env perl

use strict;
use warnings;

# Usage g_group.pl <ref_Allele_history> <from version> <to version>
# reads input from STDIN and writes converted output to STDOUT

my $ref = $ARGV[0];
my $from = $ARGV[1];
my $to = $ARGV[2];

my $ref_fh;
my %ref_in;
open($ref_fh, "<", $ref) or die "Could not open allele_history $ref";

my @header = split(/\t/,<$ref_fh>);

my $i = 0;
my %version_id;
foreach my $version (@header) {
	$version_id{$version} = $i++;
}
my $from_id = $version_id{$from};
my $to_id = $version_id{$to};
my %translate;
while (<$ref_fh>) {
	my @a = split /\t/;
	$translate{$a[$from_id]} = $a[$to_id];
}

my $input_fh;
open($input_fh, "<-") or die "Cannot open stdin";
while (<$input_fh>) {
    chomp;
    if (exists($translate{$_})) {
        print $translate{$_}."\n";
    } else {
	die "Warning unknown allele: $_\n"
    }
}
close $input_fh;

