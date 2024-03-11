library(ggplot2)
library(zoo)

## input parameters
data_dir_polyA = '/Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/'
data_dir_ori = '/Users/pmonsieurs/programming/leishmania_snsseq/results/927/'
data_dir_ori_shuffled = '/Users/pmonsieurs/programming/leishmania_snsseq/results/927/'

## get files for the polyA results
window = 2000
poly = 4
polyA_files = list.files(data_dir_polyA) #, pattern=paste0("927*.window",window,".poly_", poly, ".csv"))
polyA_files = polyA_files[grep(paste0("window",window,".poly_", poly, ".csv"), polyA_files)]
polyA_files = polyA_files[grep("927", polyA_files)]
polyA_files = polyA_files[-grep("667", polyA_files)]
polyA_files = polyA_files[-grep("668", polyA_files)]
polyA_files

## get files for the G4 hunter predictions
parameter_setting = 'merged'
cov_files_sns = list.files(data_dir_ori, pattern="*.cov")
cov_files_sns = cov_files_sns[grep(parameter_setting, cov_files_sns)]
cov_files_sns = cov_files_sns[grep("927_G4_", cov_files_sns)]

cov_files_shuffled = list.files(data_dir_ori_shuffled, pattern="*.cov")
cov_files_shuffled = cov_files_shuffled[grep("seed666", cov_files_shuffled)]
cov_files_shuffled = cov_files_shuffled[grep("927_G4_", cov_files_shuffled)]
cov_files_shuffled


## get the DRIP-seq data
drip_files =  list.files(data_dir_ori, pattern="*.cov")
drip_files = drip_files[grep("927-DRIP", drip_files)]
drip_files




# cov_files = c(cov_files_sns, cov_files_shuffled)

cov_data_all = data.frame()
first_file = 1

## first read in for the true data (G4Hunter results
## based on the SNSseq data). Afterwards for the shuffled data
for (cov_file in cov_files_sns) {
  cov_data = read.csv(paste0(data_dir_ori, cov_file), header=FALSE)
  colnames(cov_data) = c('position', 'coverage')
  cov_data = cov_data[-nrow(cov_data),]
  
  ## extract the strand from the sample name
  strand = unlist(strsplit(cov_file, split="\\."))[1]
  strand = unlist(strsplit(cov_file, split="_"))[3]
  strand = tolower(strand)
  strand = gsub("minus", 'min', strand)
  print(strand)
  cov_data$strand = strand
  cov_data$pattern = "G4"
  cov_data$pattern2 = paste0("G4 ", strand)
  
  
  
  window_size <- 100
  cov_data$coverage_smoothed = rollapply(cov_data$coverage, width = 2 * window_size + 1, FUN = mean, align = "center", fill = NA)
  
  ## create sample name by excluding the plus and min
  ## information from the strand
  sample = unlist(strsplit(cov_file, split="\\."))[2]
  sample = unlist(strsplit(sample, split="_"))[3]
  # sample = gsub("merged_", "", sample)
  sample
  
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


## read in the cov data based from G4Hunter based on the
## shuffled regions
for (cov_file in cov_files_shuffled) {
  cov_data = read.csv(paste0(data_dir_ori_shuffled, cov_file), header=FALSE)
  colnames(cov_data) = c('position', 'coverage')
  cov_data = cov_data[-nrow(cov_data),]
  
  ## extract the strand from the sample name
  strand = unlist(strsplit(cov_file, split="\\."))[1]
  strand = unlist(strsplit(cov_file, split="_"))[3]
  strand = tolower(strand)
  strand = gsub("minus", 'min', strand)
  
  print(strand)
  cov_data$strand = strand
  cov_data$pattern = "G4"
  cov_data$pattern2 = paste0("G4 ", strand)
  
  
  window_size <- 100
  cov_data$coverage_smoothed = rollapply(cov_data$coverage, width = 2 * window_size + 1, FUN = mean, align = "center", fill = NA)
  
  ## create sample name by excluding the plus and min
  ## information from the strand
  sample_full = unlist(strsplit(cov_file, split="\\."))[2]
  sample = unlist(strsplit(sample_full, split="_"))[4]
  sample = paste0("shuffled control ", sample)
  sample

  
  # sample
  cov_data$sample = sample
  
  ## add some additional column to allow visualising them
  ## together with the random / shuffled ORI
  seed = unlist(strsplit(sample_full, split="_"))[2]
  seed = gsub("seed", "", seed)
  cov_data$type = 'shuffled'
  cov_data$seed = seed
  
  if (first_file == 1) {
    cov_data_all = cov_data
    first_file = 0
  }else{
    cov_data_all = rbind.data.frame(cov_data_all, cov_data)
  }
  
}



## read in the polyA data
sample_count = 0
first_file = 1
for (nucl_file in polyA_files) {
  
  sample_count = sample_count + 1
  cov_data = read.csv(paste0(data_dir_polyA,nucl_file), header=TRUE)
  cov_data[,1] = cov_data[,1] - 2000

  
  
  ## extract the strand from the sample name, and assign to 
  ## separate column
  
  # code for normal samples
  if (grepl("seed", nucl_file)) {
    sample_data = unlist(strsplit(nucl_file, split="_"))[4]
    sample = paste0("shuffled control ", sample_data)
  }else{
    sample_data = unlist(strsplit(nucl_file, split="_"))[3]
    sample = gsub("merged_", "", sample_data)
    sample    
  }

  
  # code for shuffled controls
  sample_data = unlist(strsplit(nucl_file, split="-"))[1]
  
  
  print(paste0(nucl_file, " ---> ", sample))  
  ## calculate the percentage instead of using counts. Should only be done when
  ## the poly value is 1, so when basically looking for GC-content like things
  cov_data_perc = cov_data
  colnames(cov_data_perc)[1] = 'pos'
  # head(cov_data_perc)
  
  ## only select the polyAAA field and subsequently
  ## melt data to be able to plot
  cov_data_perc = cov_data_perc[,c('pos', 'AAAA', 'TTTT')]
  cov_data_melt = melt(cov_data_perc, id.vars = "pos")
  colnames(cov_data_melt) = c('pos', 'pattern', 'cov')
  cov_data_melt$sample = sample
  # head(cov_data_melt)

  ## do smoothing  
  window_size <- 100
  cov_data_melt$cov_smoothed = rollapply(cov_data_melt$cov, width = 2 * window_size + 1, FUN = mean, align = "center", fill = NA)
  
  ## specify the pattern, and get the strand
  cov_data_melt$strand = NA
  cov_data_melt$pattern2 = NA
  
  cov_data_melt[cov_data_melt$pattern == "AAAA",]$strand = 'plus'
  cov_data_melt[cov_data_melt$pattern == "AAAA",]$pattern2 = 'polyA plus'
  
  cov_data_melt[cov_data_melt$pattern == "TTTT",]$strand = 'min'
  cov_data_melt[cov_data_melt$pattern == "TTTT",]$pattern2 = 'polyA min'
  
  cov_data_melt$pattern = 'poly A'
  
  ## add to the overall data
  if (first_file == 1) {
    cov_data_polyA_all = cov_data_melt
    first_file = 0
  }else{
    cov_data_polyA_all = rbind.data.frame(cov_data_polyA_all, cov_data_melt)
  }
}





## read in DRIPseq data
cov_data_drip_all = data.frame()
first_file = 1

for (drip_file in drip_files) {
  cov_data = read.csv(paste0(data_dir_ori, drip_file), header=FALSE)
  colnames(cov_data) = c('position', 'coverage')
  cov_data = cov_data[-nrow(cov_data),]
  
  ## extract the strand from the sample name
  # strand = unlist(strsplit(mnase_file, split="\\."))[2]
  # strand = unlist(strsplit(strand, split="_"))[3]
  # strand = tolower(strand)
  # strand = gsub("minus", 'min', strand)
  # print(strand)
  strand = NA
  
  cov_data$strand = strand
  cov_data$pattern = "DRIP-Seq"
  cov_data$pattern2 = "DRIP-seq"
  
  
  
  window_size <- 100
  cov_data$coverage_smoothed = rollapply(cov_data$coverage, width = 2 * window_size + 1, FUN = mean, align = "center", fill = NA)
  
  ## create sample name by excluding the plus and min
  ## information from the strand
  sample = unlist(strsplit(drip_file, split="\\."))[2]
  if (grepl("seed", sample)) {
    sample = unlist(strsplit(sample, split="_"))[4]
    sample = paste0("shuffled control ", sample)
  }else{
    sample = unlist(strsplit(sample, split="_"))[3]
  }
  # sample = gsub("merged_", "", sample)
  sample
  
  # sample
  cov_data$sample = sample
  
  ## add some additional column to allow visualising them
  ## together with the random / shuffled ORI
  # cov_data$type = NA
  # cov_data$seed = NA
  
  if (first_file == 1) {
    cov_data_drip_all = cov_data
    first_file = 0
  }else{
    cov_data_drip_all = rbind.data.frame(cov_data_drip_all, cov_data)
  }
  
}



## merge both data types. First select relevant
## columns and rename column names
cov_data_all_sub = cov_data_all[,c('position', 
                                   'coverage', 
                                   'coverage_smoothed',
                                   'strand',
                                   'pattern',
                                   'pattern2',
                                   'sample')]

colnames(cov_data_all_sub) = c('pos', 'cov', 'cov_smoothed', 'strand', 'pattern', 'pattern2', 'sample')
head(cov_data_all_sub)


## do the same for the polyA information
cov_data_polyA_all_sub = cov_data_polyA_all[,c('pos',
                                               'cov',
                                               'cov_smoothed',
                                               'strand',
                                               'pattern',
                                               'pattern2',
                                               'sample')]
head(cov_data_polyA_all_sub)


cov_data_drip_all_sub = cov_data_drip_all[,c('position', 
                                               'coverage', 
                                               'coverage_smoothed',
                                               'strand',
                                               'pattern',
                                               'pattern2',
                                               'sample')]
colnames(cov_data_drip_all_sub) = c('pos', 'cov', 'cov_smoothed', 'strand', 'pattern', 'pattern2', 'sample')


## merge both data types
plot_data = rbind.data.frame(cov_data_all_sub, cov_data_polyA_all_sub)
plot_data = rbind.data.frame(plot_data, cov_data_drip_all_sub)
head(plot_data)

## create plot with separate color per line
p = ggplot(data=plot_data, aes(x=pos, y=cov_smoothed)) + 
  geom_line(aes(color=pattern2)) + 
  theme_bw() + 
  facet_wrap(~ sample) +
  ylab("") + xlab("") +
  scale_x_continuous(
    breaks = c(-1850, 0, window-150),
    labels = c("-2kb", "center", "+2kb")
  ) + 
  theme_bw() + 
  facet_wrap(~ sample) + 
  # scale_color_manual(values = colors) + 
  theme(legend.title=element_blank())

p

output_file = paste0(data_dir_polyA, parameter_setting, "_with_polyA.png")
ggsave(output_file, p, width=10, height=6)



## create plot with separate color per strand and 
## separate linetype per analysis method
p = ggplot(data=plot_data, aes(x=pos, y=cov_smoothed)) + 
  geom_line(aes(color=strand, linetype=pattern)) + 
  theme_bw() + 
  facet_wrap(~ sample) +
  ylab("") + xlab("") +
  scale_x_continuous(
    breaks = c(-1850, 0, window-150),
    labels = c("-2kb", "center", "+2kb")
  ) + 
  theme_bw() + 
  facet_wrap(~ sample) + 
  # scale_color_manual(values = colors) + 
  theme(legend.title=element_blank())

p

output_file = paste0(data_dir_polyA, parameter_setting, "_with_polyA.variant.png")
ggsave(output_file, p, width=10, height=6)






## try out with two axes. First add additional column 
## containing the polyA values. This make take some time!
head(cov_data_all_sub)
cov_data_merged = cov_data_all_sub
cov_data_merged$polyA = 0
cov_data_merged$DRIPseq = 0

for (i in 1:dim(cov_data_merged)[1]) {
  if (i%%100 == 0) {print(i)}
  sample = cov_data_merged[i,]$sample
  pos = cov_data_merged[i,]$pos
  strand = cov_data_merged[i,]$strand
  
  ## add additional polyA column to the cov_data_merged data frame
  matches = sum(cov_data_polyA_all_sub$pos == pos & cov_data_polyA_all_sub$sample == sample & cov_data_polyA_all_sub$strand == strand)
  # print(matches)
  polyA_value =  cov_data_polyA_all_sub[cov_data_polyA_all_sub$pos == pos & cov_data_polyA_all_sub$sample == sample & cov_data_polyA_all_sub$strand == strand,]$cov_smoothed
  if (is.numeric(polyA_value) && length(polyA_value) > 0) {
    # print(polyA_value)
    cov_data_merged[i,]$polyA = polyA_value
  }
  
  ## add additional DRIP-seq column to the cov_data_merged data frame
  matches = sum(cov_data_drip_all_sub$pos == pos & cov_data_drip_all_sub$sample == sample)
  # print(matches)
  drip_value =  cov_data_drip_all_sub[cov_data_drip_all_sub$pos == pos & cov_data_drip_all_sub$sample == sample,]$cov_smoothed
  if (is.numeric(mmnase_value) && length(mmnase_value) > 0) {
    # print(polyA_value)
    cov_data_merged[i,]$DRIPseq = drip_value
  }
  
  # if (i > 1000) {break}
  
}

library(scales)
cov_data_merged$DRIPseq_color = "DRIP-Seq"
manual_colors <- hue_pal()(3)
manual_colors = c(manual_colors[3], manual_colors[1:2])

## save cov_data_merged to the output diretory
# cov_data_file = paste0(data_dir_ori, 'cov_data_927.rds')
# saveRDS(cov_data, cov_data_file)

p = ggplot(data=cov_data_merged, aes(x=pos, y=cov_smoothed)) + 
  geom_line(aes(color=pattern2, y=cov_smoothed, linetype="G4"), linewidth=0.80) + 
  # geom_line(aes(x=pos,y=polyA*0.30), linewidth=0.80, color = "#FFA500") + 
  # geom_line(aes(x=pos,y=polyA*2, linetype = pattern2), linewidth=0.80, color = "#FFA500") + 
  geom_line(aes(x=pos,y=polyA, color = pattern2, linetype="polyA"), linewidth=0.80) + 
  geom_line(aes(x=pos,y=DRIPseq*0.25, color=DRIPseq_color, linetype="DRIP-Seq"), linewidth=0.80) +
  facet_wrap(~ sample) + 
  scale_y_continuous(name = "G4 and polyA", sec.axis = sec_axis(~./0.25, name = "DRIP-Seq")) + 
  xlab("") +
  scale_x_continuous(
    breaks = c(-1850, 0, window-150),
    labels = c("-2kb", "center", "+2kb")) +
  theme_bw() + 
  scale_color_manual(values = manual_colors) + 
  theme(panel.spacing = unit(0.5, "cm"),
        legend.title=element_blank(),
        legend.key.width = unit(1.5, "cm")) + 
  # guides(linetype = guide_legend(override.aes = list(linetype = c("G4" = "solid", "polyA" = "dashed"))))
  scale_linetype_manual(values = c("G4" = "solid", "polyA" = "dashed", "DRIP-Seq" = "solid"))

p



output_file = paste0(data_dir_ori, "927_with_polyA.variant.dual_axes.tiff")
ggsave(output_file, p, width=10, height=6, dpi=300)


## find the maximum hight of the peak of the MNA-seq data verus 
## the ORI position (= position 0), this mean finding the maximum
## value for the DRIPseq on the minus region, and on the plus region


for (sample in unique(cov_data_merged$sample)) {
  
  print(paste0("sample --> ", sample))
  cov_data_max = cov_data_merged[cov_data_merged$sample == sample,]
  
  ## focus on downstream side
  print("..downstream")
  cov_data_max_down = cov_data_max[cov_data_max$pos < 0,]
  max_DRIPseq_index = which(cov_data_max_down$DRIPseq == max(cov_data_max_down$DRIPseq, na.rm = TRUE))
  max_DRIPseq_peak = cov_data_max_down[max_DRIPseq_index,]$pos[1]
  print(paste0("....DRIPseq peak: ", max_DRIPseq_peak))
  
  ## split between plus and minus strand
  print("......strand +")
  cov_data_max_down_plus = cov_data_max_down[cov_data_max_down$strand == "plus",] 
  
  max_polyA_index = which(cov_data_max_down_plus$polyA == max(cov_data_max_down_plus$polyA, na.rm = TRUE))
  max_polyA_peak = cov_data_max_down_plus[max_polyA_index,]$pos[1]
  print(paste0("........polyA peak: ", max_polyA_peak))
  
  max_G4_index = which(cov_data_max_down_plus$cov_smoothed == max(cov_data_max_down_plus$cov_smoothed, na.rm = TRUE))
  max_G4_peak = cov_data_max_down_plus[max_G4_index,]$pos[1]
  print(paste0("........G4 peak: ", max_G4_peak))
  
  print("......strand -")
  cov_data_max_down_min = cov_data_max_down[cov_data_max_down$strand == "min",] 
  max_polyA_index = which(cov_data_max_down_min$polyA == max(cov_data_max_down_min$polyA, na.rm = TRUE))
  max_polyA_peak = cov_data_max_down_min[max_polyA_index,]$pos[1]
  print(paste0("........polyA peak: ", max_polyA_peak))
  
  max_G4_index = which(cov_data_max_down_min$cov_smoothed == max(cov_data_max_down_min$cov_smoothed, na.rm = TRUE))
  max_G4_peak = cov_data_max_down_min[max_G4_index,]$pos[1]
  print(paste0("........G4 peak: ", max_G4_peak))
  
  ## focus on upstream side
  print("...upstream")
  cov_data_max_up = cov_data_max[cov_data_max$pos > 0,]
  max_DRIPseq_index = which(cov_data_max_up$DRIPseq == max(cov_data_max_up$DRIPseq, na.rm = TRUE))
  max_DRIPseq_peak = cov_data_max_up[max_DRIPseq_index,]$pos[1]
  print(paste0("...DRIPseq peak: ", max_DRIPseq_peak))
  
  
  ## split between plus and minus strand
  print("......strand +")
  cov_data_max_up_plus = cov_data_max_up[cov_data_max_up$strand == "plus",] 
  
  max_polyA_index = which(cov_data_max_up_plus$polyA == max(cov_data_max_up_plus$polyA, na.rm = TRUE))
  max_polyA_peak = cov_data_max_up_plus[max_polyA_index,]$pos[1]
  print(paste0("........polyA peak: ", max_polyA_peak))
  
  max_G4_index = which(cov_data_max_up_plus$cov_smoothed == max(cov_data_max_up_plus$cov_smoothed, na.rm = TRUE))
  max_G4_peak = cov_data_max_up_plus[max_G4_index,]$pos[1]
  print(paste0("........G4 peak: ", max_G4_peak))
  
  print("......strand -")
  cov_data_max_up_min = cov_data_max_up[cov_data_max_up$strand == "min",] 
  max_polyA_index = which(cov_data_max_up_min$polyA == max(cov_data_max_up_min$polyA, na.rm = TRUE))
  max_polyA_peak = cov_data_max_up_min[max_polyA_index,]$pos[1]
  print(paste0("........polyA peak: ", max_polyA_peak))
  
  max_G4_index = which(cov_data_max_up_min$cov_smoothed == max(cov_data_max_up_min$cov_smoothed, na.rm = TRUE))
  max_G4_peak = cov_data_max_up_min[max_G4_index,]$pos[1]
  print(paste0("........G4 peak: ", max_G4_peak))
  
  
  
  
}





# Load the ggplot2 package
library(ggplot2)

# Create a sample dataframe
data <- data.frame(
  x = 1:10,
  y1 = c(3, 6, 8, 12, 10, 9, 7, 5, 4, 2),
  y2 = c(50, 45, 40, 35, 30, 25, 20, 15, 10, 5)
)

# Create a ggplot with two different y-axes
ggplot(data, aes(x = x)) +
  geom_line(aes(y = y1), color = "blue") +
  geom_line(aes(y = y2 * 0.1), color = "red") +
  scale_y_continuous(name = "Primary Y-Axis", sec.axis = sec_axis(~.*10, name = "Secondary Y-Axis (scaled by 10)")) +
  labs(title = "GGPlot with Two Different Y-Axes")




