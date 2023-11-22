library(ggplot2)
library(zoo)
library(reshape)

## get all the coverage files for the normal ORI files
data_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/seqkit/'
setwd(data_dir)


## type is the number of nucleotides
type = 'mono'
type = 'tetra'
nucl_files = list.files(data_dir, pattern=paste0("*", type, ".txt"))

nucl_data_all = data.frame()
first_file = 1



for (nucl_file in nucl_files) {
  
  nucl_data = read.table(nucl_file, header=FALSE)
  
  
  ## extract the strand from the sample name, and assign to 
  ## separate column
  sample_data = unlist(strsplit(nucl_file, split="_"))[1:3]
  sample_name= paste0(sample_data[1], "_", sample_data[2], "_", sample_data[3])
  print(sample_name)
  nucl_data$sample = sample_name
  
  ## do for specific 
  
  ## only select the relevant columns
  nucl_data = nucl_data[,c(2,3,4,5,6,8)]
  colnames(nucl_data) = c('pattern_name', 'pattern', 'strand', 'start', 'end',  'sample')  
  
  ## get the coverage of a certain pattern per position
  ## using the "table" function, in this case using the 
  ## xtabs function selected on two columns
  cov_data = t(xtabs(~ pattern + start, nucl_data))
  head(cov_data)
  
  ## calculate the percentage instead of using counts
  cov_data_perc = cov_data
  cov_data_perc = cov_data_perc[,1:4]/rowSums(cov_data_perc[,1:4])
  
  ## melt data to be able to plot
  cov_data_melt = melt(cov_data_perc)
  colnames(cov_data_melt) = c('pos', 'pattern', 'cov')
  cov_data_melt$sample = sample_name
  head(cov_data_melt)
  

  
  window_size <- 100
  cov_data_melt$cov_smoothed = rollapply(cov_data_melt$cov, width = 2 * window_size + 1, FUN = mean, align = "center", fill = NA)
  
  ## specify the pattern
  # nucl_data = nucl
  
  ## add to the overall data
  if (first_file == 1) {
    cov_data_all = cov_data_melt
    first_file = 0
  }else{
    cov_data_all = rbind.data.frame(cov_data_all, cov_data_melt)
  }
}


# cov_data_all = cov_data_all[cov_data_all$sample == 'merged_BSF-b_ORIs',]


ggplot(data=cov_data_all, aes(x=pos, y=cov_smoothed)) + 
  geom_line(aes(color=pattern)) + 
  coord_cartesian(xlim=c(300,4700), ylim=c(0, 0.50)) +
  theme_bw() + 
  facet_wrap(~ sample, scales = "free_y") 
