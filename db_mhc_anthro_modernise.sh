#!/bin/bash

set -euo pipefail

basedir=`dirname $0`

inputfile=$1

(head -n1 ${inputfile};
paste \
<(cut -f 1-10 ${inputfile} | tail -n+2) \
<(cut -f 11 ${inputfile} | tail -n+2 | sed 's/^/HLA-A*/' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1050 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt ) \
<(cut -f 12 ${inputfile} | tail -n+2 | sed 's/^/HLA-A*/' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1050 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 13 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-Cw*\1/' -e 's/HLA-Cw\*0701/HLA-Cw*07011/g' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1080 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 14 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-Cw*\1/' -e 's/HLA-Cw\*0701/HLA-Cw*07011/g' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1080 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 15 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-B*\1/' -e 's/HLA-B\*4701/HLA-B*4701101/g' -e 's/HLA-B\*1517/HLA-B*1517101/g' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1110 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 16 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-B*\1/' -e 's/HLA-B\*4701/HLA-B*4701101/g' -e 's/HLA-B\*1517/HLA-B*1517101/g' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1110 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 17-20 ${inputfile} | tail -n+2) \
<(cut -f 21 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-DQA1*\1/' -e 's/^HLA-DQA1\*0102/HLA-DQA1*01021/g' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1090 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 22 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-DQA1*\1/' -e 's/^HLA-DQA1\*0102/HLA-DQA1*01021/g' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1090 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 23 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-DQB1*\1/' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1090 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 24 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-DQB1*\1/' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1090 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 25 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-DPA1*\1/' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1090 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 26 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-DPA1*\1/' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1090 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 27 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-DPB1*\1/' -e 's/^HLA-DPB1\*3400/HLA-DPB1*3401/g' -e 's/^HLA-DPB1\*0400/HLA-DPB1*0401/g' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1090 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
<(cut -f 28 ${inputfile} | tail -n+2 | sed -e 's/^\([0-9]\+\)/HLA-DPB1*\1/' -e 's/^HLA-DPB1\*3400/HLA-DPB1*3401/g' -e 's/^HLA-DPB1\*0400/HLA-DPB1*0401/g' | ${basedir}/hla_version_convert.pl Allelelist_history.txt 1090 3630 1 0 | ${basedir}/g_group.pl hla_nom_g.txt) \
)
