library(ggplot2)

data_dir = '/Users/pmonsieurs/programming/trypanosoma_sofi_mining/results/ori/'
setwd(data_dir)

cov_files = list.files(data_dir, pattern="*.cov")

cov_data_all = data.frame()
first_file = 1

for (cov_file in cov_files) {
  cov_data = read.csv(cov_file, header=FALSE)
  colnames(cov_data) = c('position', 'coverage')
  cov_data = cov_data[-nrow(cov_data),]
  cov_data$sample = cov_file
  
  if (first_file == 1) {
    cov_data_all = cov_data
    first_file = 0
  }else{
    cov_data_all = rbind.data.frame(cov_data_all, cov_data)
  }
  
}

ggplot(data=cov_data_all, aes(x=position, y=coverage)) + 
  geom_line(aes(group=sample, colour=sample)) + 
  theme_bw()

ggplot(data=cov_data_all, aes(x=position, y=coverage)) + 
  geom_line() + 
  theme_bw() + 
  facet_wrap(~ sample)


