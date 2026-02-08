#!/usr/bin/env perl

#  hla_summary.pl
#  hla_db_util - set of utilities for working with IMGT/HLA database data
#
# Created by Martin Pollard.
# Copyright © 2017, 2026 Genome Research Limited.
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

use HLA qw(print_full_matches parse_flat print_full_matches_multi);

my $ref = $ARGV[0];
my $input = $ARGV[1];
my $output_full = $ARGV[2];
my $output_one = $ARGV[3];
my $error_name = $ARGV[4];

my %validation_data_hash = parse_flat($ref);
my %input_data = parse_flat($input);

print_full_matches($output_full,$error_name,%validation_data_hash, %input_data);
