library(ggplot2)

src_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/'
data_dir = paste0(src_dir, 'data/')

plot_data = data.frame()
for (i in 1:4) {
  file_name = paste0('read_length_count.', i, '_S', i, '.csv')
  file_name = paste0(data_dir, file_name)
  data = read.csv(file_name, sep=" ", header=FALSE)
  colnames(data) = c('readL', 'count')
  data$sample = paste0('S', i)
  plot_data = rbind(plot_data, data)
}


p = ggplot(plot_data, aes(x=readL, y=log10(count)))
p = p + geom_line(aes(colour=sample))
p


# calculate relative percentages
for (i in 1:4) {
  sample = paste0('S', i)
  print(sample)
  short_count = sum(plot_data[plot_data$sample == sample & plot_data$readL < 100,]$count)
  percentage_short = short_count/sum(plot_data[plot_data$sample == sample,]$count)
  print(percentage_short)
}

