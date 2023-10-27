library(ggplot2)
library(zoo)
library(reshape)

## get all the coverage files for the normal ORI files
data_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/'
setwd(data_dir)


## type is the number of nucleotides
# type = 'mono'
# type = 'tetra'

poly = 8
window=2000
nucl_files = list.files(data_dir, pattern=paste0("*.window",window,".poly_", poly, ".csv"))
nucl_files = nucl_files[-grep("667", nucl_files)]
nucl_files = nucl_files[-grep("668", nucl_files)]

nucl_data_all = data.frame()
first_file = 1

sample_names = c("BSF ORIs", "PCF ORIs", "shuffled control", "shuffled control 2")
sample_count = 0

for (nucl_file in nucl_files) {
  
  sample_count = sample_count + 1
  cov_data = read.csv(nucl_file, header=TRUE)
  cov_data[,1] = cov_data[,1] + 1
  
  
  ## extract the strand from the sample name, and assign to 
  ## separate column
  sample_data = unlist(strsplit(nucl_file, split="_"))[1:3]
  sample_name= paste0(sample_data[1], "_", sample_data[2], "_", sample_data[3])
  print(sample_name)
  
  sample_name = sample_names[sample_count]
  
  

  ## calculate the percentage instead of using counts. Should only be done when
  ## the poly value is 1, so when basically looking for GC-content like things
  cov_data_perc = cov_data
  if (poly == 1) {
    cov_data_perc[,2:5] = cov_data_perc[,2:5]/rowSums(cov_data_perc[,2:5])
  }
  colnames(cov_data_perc)[1] = 'pos'
  head(cov_data_perc)
  
  ## melt data to be able to plot
  cov_data_melt = melt(cov_data_perc, id.vars = "pos")
  head(cov_data_melt)
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


## define the colors
colors <- c("darkgreen", "darkred", "#FFA500", "darkblue")

p = ggplot(data=cov_data_all, aes(x=pos, y=cov_smoothed)) + 
  geom_line(aes(color=pattern), linewidth=0.80) + 
  coord_cartesian(xlim=c(250,2*window-250)) +#, ylim=c(0, 0.75)) +
  ylab("") + xlab("") +
  scale_x_continuous(
    breaks = c(150, window, 2*window-150),
    labels = c("-2kb", "center", "+2kb")
  ) + 
  theme_bw() + 
  facet_wrap(~ sample) + 
  scale_color_manual(values = colors) + 
  theme(panel.spacing = unit(0.5, "cm"),
        legend.title=element_blank())

p

output_file = paste0(data_dir, 'poly_nucleotide.poly', poly, '.png')
ggsave(file = output_file, plot=p, width=10, height=6)

output_file = paste0(data_dir, 'poly_nucleotide.poly', poly, '.tiff')
# ggsave(file = output_file, plot=p, dpi=300, width=9, height=6)


