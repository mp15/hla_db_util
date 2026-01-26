#!/usr/bin/env Rscript

#  align.R
#  hla_db_util - set of utilities for working with IMGT/HLA database data
#
# Created by Martin Pollard.
# Copyright © 2018 Genome Research Limited.
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

# Plot variation across alignments

library('seqinr')
library(ggplot2)

plot.variationvec <- function(gene){
    align <-read.alignment(paste0(gene,'_gen.msf'),'msf')
    align_matrix <- as.matrix(align)
    #dist_align <- dist.alignment(align)

    allele.vector <- function(align.matrix) {
        apply(align_matrix,2, function(vec) {length(table(vec[vec!='-']))-1})
    }
    variation_vec <- allele.vector(align_matrix)
    #For some reason there are sites where we have no data in the alignment?
    variation_vec[variation_vec==-1] <- NA
    vp <- ggplot(data.frame(Position=seq_along(variation_vec),Variants=variation_vec), aes(Position,Variants)) + geom_bar(stat="identity") + ggtitle(toupper(gene)) + coord_fixed(ratio = 500)

    ggsave(paste0(gene,'_variations.pdf'),vp, height=4)
}
plot.variationvec_alleles <- function(gene){
    align <-read.alignment(paste0(gene,'_gen.msf'),'msf')
    alleles <- unlist(read.table('A.allele_list'))
    align_matrix <- as.matrix(align)
    align_matrix <- align_matrix[alleles,]
    #dist_align <- dist.alignment(align)

    allele.vector <- function(align.matrix) {
        apply(align_matrix,2, function(vec) {length(table(vec[vec!='-']))-1})
    }
    variation_vec <- allele.vector(align_matrix)
    #For some reason there are sites where we have no data in the alignment?
    variation_vec[variation_vec==-1] <- NA
    vp <- ggplot(data.frame(Position=seq_along(variation_vec),Variants=variation_vec), aes(Position,Variants)) + geom_bar(stat="identity") + ggtitle(toupper(gene) + coord_fixed(ratio = 500))
    ggsave(paste0(gene,'_variations_gdap.pdf'),vp, height=4)
}

#plot number of vario
for (i in c('a','b','c','dpa1','dpb1','dqa1','dqb1','drb1', 'drb3','drb4','drb5')) {
    plot.variationvec(i)
}

# Perhaps plot similarly % as heatmap?

#library(RColorBrewer)
#library(reshape2)
#ggplot(dist_align) +   geom_tile()
#image(as.matrix(dist_align), xlab = 'Matrix rows', ylab = 'Matrix columns',axes = F))

#dist_align.melted <- melt(as.matrix(dist_align))
#ggplot(dist_align.melted, aes(x = Var1, y = Var2, fill = value)) + geom_tile()

for (i in c('a','b','c','dpb1','dqa1','dqb1','drb1')) {
    plot.variationvec_alleles(i)
}
