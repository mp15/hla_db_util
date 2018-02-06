library('seqinr')
library(ggplot2)

plot.variationvec <- function(gene){
    align <-read.alignment(paste0(gene,'_gen.msf'),'msf')
    align_matrix <- as.matrix(align)
    dist_align <- dist.alignment(align)

    allele.vector <- function(align.matrix) {
        apply(align_matrix,2, function(vec) {length(table(vec[vec!='-']))-1})
    }
    variation_vec <- allele.vector(align_matrix)
    #For some reason there are sites where we have no data in the alignment?
    variation_vec[variation_vec==-1] <- NA
    vp <- ggplot(data.frame(Position=seq_along(variation_vec),Variants=variation_vec), aes(Position,Variants)) + geom_bar(stat="identity") + ggtitle(toupper(gene))
    ggsave(paste0(gene,'_variations.pdf'),vp)
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
