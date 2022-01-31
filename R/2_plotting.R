library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)

date = Sys.Date()
load('data/sc.Robj')

#load metadata
url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FFinalAnnotationTable%2Etxt%2Egz'
tmp <- tempfile()
##
download.file(url,tmp)
metadata <- read.delim(gzfile(tmp), sep = '\t', stringsAsFactors = F)
colnames(metadata)[1] <- 'ID'

#build metadata data frame
df <- data.frame('Cluster' = sc@cpart, sc@tsne) 
df$ID <- gsub('\\.', ':', rownames(df))

df <- df %>% 
  inner_join(metadata)

df$Condition <- factor(df$Condition, levels = c('Ctrl', 'MS'))

df$Celltypes <- factor(df$Celltypes, levels = c('OPCs', 'COPs', 'ImOlGs', 'Oligo1', 'Oligo2', 'Oligo3', 'Oligo4', 'Oligo5', 'Oligo6', 'Astrocytes', 'Astrocytes2', 'Neuron1','Neuron2','Neuron3', 'Neuron4', 'Neuron5', 'Endothelial_cells1', 'Endothelial_cells2', 'Macrophages', 'Microglia_Macrophages', 'Immune_cells', 'Vasc_smooth_muscle','Pericytes'))

df$Lesion <- factor(df$Lesion, levels = c('Ctrl', 'NAWM', 'A', 'CA', 'CI', 'RM'))

cell_numbers <- numeric()
for (i in 1:max(sc@cpart, na.rm = T)) {
  cell_numbers[i] <- length(na.omit(sc@cpart[sc@cpart==i]))
}
names(cell_numbers) <- c(1:max(sc@cpart, na.rm = T))
retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(sc@ndata)[2]/100]))

df <- df[df$Cluster %in% retain_cl,]

#export means per celltype gene expression table
colnames(sc@ndata) <- gsub('\\.', ':', colnames(sc@ndata))

#plot igfr expression
#violin plot of gene expressions
igf1r <- data.frame('IGF1R' = as.matrix(sc@ndata)['IGF1R', ]* min(sc@counts))
rownames(igf1r) <- gsub('\\.', ':', rownames(igf1r))
data <- cbind(df, igf1r[df$ID,])
colnames(data)[11] <- 'IGF1R' 


sorted <- data %>% group_by(Lesion, Celltypes) %>%
  summarise(expression = mean(IGF1R)) %>%
  arrange(expression)
sorted2 <- data %>% group_by(Celltypes) %>%
  summarise(expression = mean(IGF1R)) %>%
  arrange(expression)


data$Celltypes <- factor(data$Celltypes, levels = as.character(sorted2$Celltypes))


ggplot(data[!data$Celltypes %in% c('Macrophages', 'Vasc_smooth_muscle'),], aes(Celltypes, IGF1R, fill=Celltypes)) +
  stat_summary(fun = mean, geom = "bar", color = 'black', lwd=0.25) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",width = 0) +
  coord_flip() +
  theme_minimal() +
  labs(y='IGF1R expression') +
  facet_wrap(~Lesion, scales = 'free_x')#+
#scale_fill_manual(values = rev(colors_many))

ggsave(paste0('plots/others/',date, '-celltype-and-lesion-dependent-igf1r-expression.pdf'))


