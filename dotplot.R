library(ggplot2)
library(reshape)

raw_data <- read.csv("~/Podg_lab/IB_fasciola/tpms.csv" , sep = "\t")
colnames(raw_data) <- c("Transcript_name", "Egg_1", "Egg_2", "Juvenile_1", "Juvenile_2", "Metacerсarea_1", "Metacerсarea_2", "Adult_1", "Adult_2")
raw_data$Transcript_name <- as.character(raw_data$Transcript_name)
raw_data$type <- c(rep("repeat", 18), rep("gene", 7))

data <- melt(raw_data)
data$replicate <- 1
data[grepl(".*_2", data$variable) , 5] <- 2
names(data)[3] <- "Stage"
names(data)[4] <- "tpm"
data$Stage <- sub("_.", "", data$Stage)

ggplot(data, aes(type, tpm)) +
  geom_boxplot(alpha = 0.7) +
  geom_dotplot(binaxis = "y", stackdir= "center", position = "dodge", aes(color = replicate, fill = type), alpha = 0.6, dotsize = 1.1) +
  
  scale_y_log10() +
  
  facet_wrap(. ~Stage, ) +
  theme_grey() %+replace% 
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA),
        strip.text.x = element_text(size = 40 / .pt), 
        axis.text.x = element_text(size = 30 / .pt),
        axis.text.y = element_text(size = 30 / .pt),
        complete = TRUE)
  
  


