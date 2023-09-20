library(ggplot2)


## package if you want to do some smoothing on the values: 
if (!requireNamespace("zoo", quietly = TRUE)) {
  install.packages("zoo")
}
library(zoo)

## get all the coverage files for the normal ORI files
data_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/ori/'
setwd(data_dir)
cov_files = list.files(data_dir, pattern="*.cov")

cov_data_all = data.frame()
first_file = 1

for (cov_file in cov_files) {
  cov_data = read.csv(cov_file, header=FALSE)
  colnames(cov_data) = c('position', 'coverage')
  cov_data = cov_data[-nrow(cov_data),]
  
  ## extract the strand from the sample name
  strand = unlist(strsplit(cov_file, split="\\."))[3]
  cov_data$strand = strand
  
  window_size <- 100
  cov_data$coverage_smoothed = rollapply(cov_data$coverage, width = 2 * window_size + 1, FUN = mean, align = "center", fill = NA)
    
  ## create sample name by excluding the plus and min
  ## information from the strand
  sample = paste(unlist(strsplit(cov_file, split="\\."))[1], 
                unlist(strsplit(cov_file, split="\\."))[2],
                unlist(strsplit(cov_file, split="\\."))[4],
                sep = ".")
  # sample
  cov_data$sample = sample
  
  ## add some additional column to allow visualising them
  ## together with the random / shuffled ORI
  cov_data$type = 'ori'
  cov_data$seed = 'ori'
  
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
  geom_line(aes(color=strand)) + 
  theme_bw() + 
  facet_wrap(~ sample)

## with smoothed values
ggplot(data=cov_data_all, aes(x=position, y=coverage_smoothed)) + 
  geom_line(aes(color=strand)) + 
  theme_bw() + 
  facet_wrap(~ sample, scales = "free_y")




## read in the data from the random controls - shuffled ORI 
## sequences produced by Bridlin

data_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/ori_shuffled/'
setwd(data_dir)
cov_files = list.files(data_dir, pattern="*.cov")

cov_data_all_shuffled = data.frame()
first_file = 1

for (cov_file in cov_files) {
  cov_data = read.csv(cov_file, header=FALSE)
  colnames(cov_data) = c('position', 'coverage')
  cov_data = cov_data[-nrow(cov_data),]
  
  # sample_name_g4 = gsub("^(.*?)\\..*?\\.", "\\1", cov_file)
  # sample_name_ori = gsub(".*_([A-Z]+)_.*", "\\1", cov_file)

  ## split on the "_" and reassemble
  parts <- unlist(strsplit(cov_file, "_"))
  parts_4 = unlist(strsplit(parts[4], "\\."))[1]
  strand = unlist(strsplit(parts[4], "\\."))[2]
  sample_name_g4 = parts[1:4]
  sample_name_ori = parts[6]
  sample_name = paste0(parts[1], "_", parts[2], "_", parts[3], "_", parts_4, ".merged_", parts[6])
    
  # print(sample_name_g4)
  # print(sample_name_ori)
  print(sample_name)    

  cov_data$sample = sample_name
  cov_data$type = 'shuffled'
  cov_data$seed = parts[5]
  cov_data$strand = strand
  
  if (first_file == 1) {
    cov_data_all_shuffled = cov_data
    first_file = 0
  }else{
    cov_data_all_shuffled = rbind.data.frame(cov_data_all_shuffled, cov_data)
  }
  
}


ggplot(data=cov_data_all_shuffled, aes(x=position, y=coverage, group=seed)) + 
  geom_line(linewidth = 0.1) + 
  theme_bw() + 
  facet_wrap(~ sample)


## merged both original and shuffled
cov_data_merged = rbind(cov_data_all, 
                        cov_data_all_shuffled)

ggplot(data=cov_data_merged, aes(x=position, y=coverage, group=interaction(seed,strand))) + 
  geom_line(aes(color = seed, linetype=strand ), linewidth = 0.3) + 
  theme_bw() + 
  facet_wrap(~ sample)

ggplot(data=cov_data_merged, aes(x=position, y=coverage, group=interaction(seed,strand))) + 
  geom_line(aes(color = type, linetype=strand), linewidth = 0.3) + 
  theme_bw() + 
  facet_wrap(~ sample)
