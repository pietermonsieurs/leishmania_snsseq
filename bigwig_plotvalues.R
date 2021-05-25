library('ggplot2')

test_wig_file = '/Users/pmonsieurs/programming/leishmania_snsseq/software/test.wig'

wig_data = read.csv(test_wig_file, sep="\t", head=FALSE, comment.char="#")
colnames(wig_data) = c('chrom', 'start', 'end', 'value')

wig_data_core = wig_data[grep("core", wig_data$chrom),]

p = ggplot(wig_data_core, aes(x=value))
p = p + geom_density(aes(colour=chrom))
p = p + xlim(-5,5)
p
