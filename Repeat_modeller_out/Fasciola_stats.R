library(data.table)
library(ggplot2)
setwd("~/R_stuff/Transposons/")
if (!file.exists("RM_rawout.zip")){
  unzip("f.hep_RM.out.zip")
  RM_rawout <- fread("fasciola_hepatica.PRJEB6687.WBPS9.genomic.fa.out" , sep = " ", blank.lines.skip = T, fill = T)
  RM_rawout <- RM_rawout[-1, -16]
  names(RM_rawout) <- c("SW_score", "div.%", "del.%", "ins.%", "scaffold_name", "query_start", "query_stop", 
                      "(left)", "strand", "element_name", "class/family", "repeat_start", "repeat_end", "(left)", "ID")
  #RM_rawout$strand <- NULL
  RM_rawout$`(left)` <- NULL
  RM_rawout$`(left)` <- NULL
  RM_rawout$query_start <- as.numeric(RM_rawout$query_start)
  RM_rawout$query_stop <- as.numeric(RM_rawout$query_stop)
  RM_rawout$`class/family` <- sapply(RM_rawout$`class/family`, function(x) strsplit(x, "/")[[1]])
  RM_rawout$class <- lapply(RM_rawout$`class/family`, function(x) x[[1]])
  RM_rawout$family <- lapply(RM_rawout$`class/family`, function(x) ifelse(length(x) > 1, x[[2]], NA))
  RM_rawout$`class/family` <- NULL
  fwrite(RM_rawout, file = "RM_rawout.tsv", sep = "\t")
  zip("RM_rawout.zip", "RM_rawout.tsv")
  file.remove("RM_rawout.tsv")
  file.remove("fasciola_hepatica.PRJEB6687.WBPS9.genomic.fa.out")
}else{
  unzip("RM_rawout.zip")
  RM_rawout <- fread("RM_rawout.tsv")
  file.remove("RM_rawout.tsv")
}

#Отфильтруем все то, что нам не нужно
RM_rawout <- RM_rawout[ class != "Simple_repeat"]
RM_rawout <- RM_rawout[ class != 'Low_complexity']
RM_rawout <- RM_rawout[ class != 'tRNA']
RM_rawout <- RM_rawout[ family != 'tRNA']
RM_rawout <- RM_rawout[ class != 'rRNA']
RM_rawout <- RM_rawout[ class != 'snRNA']
RM_rawout <- RM_rawout[ class != 'srpRNA']
RM_rawout <- RM_rawout[ class != 'scRNA']
RM_rawout <- RM_rawout[ class != 'ARTEFACT']
RM_rawout$class <- unlist(RM_rawout$class)
RM_rawout$family <- unlist(RM_rawout$family)
RM_rawout$class <- sapply(RM_rawout$class, function(x) sub("?", replacement = "", x, fixed = T))
RM_rawout$family <- sapply(RM_rawout$family,  function(x) sub("?", replacement = "", x, fixed = T))
RM_rawout$element_name <- sapply(RM_rawout$element_name,  function(x) sub("?", replacement = "", x, fixed = T))
unique(RM_rawout$family)

#Делаем барплот по распределению
#Для начала - самые распространненные элементы по именам

by_elements <- as.data.table(table(RM_rawout$element_name))
by_class <- as.data.table(table(RM_rawout$class))
by_family <- as.data.table(table(RM_rawout$family))

by_elements <- by_elements[order(by_elements$N, decreasing = T), ]
by_class <- by_class[order(by_class$N, decreasing = T), ]
by_family <- by_family[order(by_family$N, decreasing = T), ]

#Plot of the top TE families

ggplot(by_family[1:20,], aes(x=reorder(V1,-N), N))+
  geom_bar(stat = "identity", colour = "black", fill = "red", alpha = "0.6")+
  labs(title = "Top20 TE families",
       x = "TE family name",
       y = "Number in genome")+
theme_grey() %+replace% 
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5),
        complete = TRUE)

#Plot of the TE classes
ggplot(by_class, aes(x=reorder(V1,-N), N))+
  geom_bar(stat = "identity", colour = "black", fill = "red", alpha = "0.6")+
  labs(title = "TE classes distribution",
       x = "TE class name",
       y = "Number in genome")+
  theme_grey() %+replace% 
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5),
        complete = TRUE)

#Plot of the elements names
ggplot(by_elements[1:15], aes(x=reorder(V1,-N), N))+
  geom_bar(stat = "identity", colour = "black", fill = "red", alpha = "0.6")+
  labs(title = "Top15 TE's",
       x = "TE name",
       y = "Number in genome")+
  theme_grey() %+replace% 
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5),
        complete = TRUE)

#Теперь графики по занимаемому проценту генома
fasciola_size=1275107963
bp_class <- RM_rawout[, .(sum(abs(query_stop - query_start))), by = class]
bp_class <- bp_class[order(bp_class$V1, decreasing = T), ]

bp_family <- RM_rawout[, .(sum(abs(query_stop - query_start))), by = family]
bp_family <- bp_family[order(bp_family$V1, decreasing = T), ]

bp_element <- RM_rawout[, .(sum(abs(query_stop - query_start))), by = element_name]
bp_element <- bp_element[order(bp_element$V1, decreasing = T), ]

#По занятым нуклеотидам
ggplot(bp_family[1:15,], aes(x=reorder(family,-V1), V1))+
  geom_bar(stat = "identity", colour = "black", fill = "green", alpha = "0.6")+
  labs(title = "Top15 TE families",
       x = "TE family name",
       y = "Occupied nucleotides")+
  theme_grey() %+replace% 
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5),
        complete = TRUE)

ggplot(bp_class, aes(x=reorder(class,-V1), V1))+
  geom_bar(stat = "identity", colour = "black", fill = "green", alpha = "0.6")+
  labs(title = "TE classes by occupancy",
       x = "TE class name",
       y = "Occupies nucleotides")+
  theme_grey() %+replace% 
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5),
        complete = TRUE)

ggplot(bp_element[1:15], aes(x=reorder(element_name,-V1), V1))+
  geom_bar(stat = "identity", colour = "black", fill = "green", alpha = "0.6")+
  labs(title = "Top15 TE elements",
       x = "TE element name",
       y = "Total nucleotides occupied")+
  theme_grey() %+replace% 
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5),
        complete = TRUE)


#Теперь красивые таблички будут с циферками

#По классам
class_table <- by_class
names(class_table) <- c("class", "Occurences")
setkey(class_table, class)
class_table <- class_table[bp_class]
names(class_table) <- c("class", "Occurences", "Total_bp")
class_table$mean_length <- class_table[, .(Mean_length = round(Total_bp/Occurences, digits = 0))]
class_table$percent_genome <- class_table$Total_bp/fasciola_size*100
fwrite(class_table, file = "graphs/classes.tsv", sep = "\t")

#По семействам
family_table <- by_family
names(family_table) <- c("family", "Occurences")
setkey(family_table, family)
family_table <- family_table[bp_family]
names(family_table) <- c("family", "Occurences", "Total_bp")
family_table$mean_length <- family_table[, .(Mean_length = round(Total_bp/Occurences, digits = 0))]
family_table$percent_genome <- family_table$Total_bp/fasciola_size*100
fwrite(family_table, file = "graphs/families.tsv", sep = "\t")

#По TE
elements_table <- by_elements
names(elements_table) <- c("elements", "Occurences")
setkey(elements_table, elements)
elements_table <- elements_table[bp_element]
names(elements_table) <- c("elements", "Occurences", "Total_bp")
elements_table$mean_length <- elements_table[, .(Mean_length = round(Total_bp/Occurences, digits = 0))]
elements_table$percent_genome <- elements_table$Total_bp/fasciola_size*100
fwrite(elements_table, file = "graphs/elements.tsv", sep = "\t")
