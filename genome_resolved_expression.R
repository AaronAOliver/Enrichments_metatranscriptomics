## Set up environment
library(ggplot2)
library(ggpmisc)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(plyr)
library(ggrepel)
library(cowplot)
library(forcats)
library(paletteer)
options(scipen = 999)

setwd("Enrichments_metatranscriptomics")

## Set graphical parameters for different ggplots
# Color mapping
taxonomy_coloring = c(
  "Alphaproteobacteria" = "#006F6C", 
  "Bacillota" = "#BC0000", 
  "Bacteroidota" = "#006DDB", 
  "Candidatus Tefagiota" = "#db6d00",
  "Desulfobacterota" = "#22CF22", 
  "Fusobacteriota" = "#A7EAF7", 
  "Gammaproteobacteria" = "#8453b5", 
  "Spirochaetota" = "#ffdf4d", 
  "Verrucomicrobiota" = "#d5b2f7"
)

# Stacked bar graphs
bar_graph_params = list(geom_bar(position="stack", stat="identity"), theme_bw(),
                        scale_fill_manual(name = "MAG Taxonomy", values = taxonomy_coloring),
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")), 
                        xlab("Metatranscriptomic Sample"))

# Beeswarm graphs
beeswarm_params = list(geom_beeswarm(cex = 0.5, size = 2), theme_bw(),
                       scale_color_manual(name = "MAG Taxonomy", values = taxonomy_coloring),
                       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")),
                       xlab("Metatranscriptomic Sample"), scale_x_discrete(drop = FALSE))

area_plot_params = list(geom_area(stat="identity", position = "stack"),
                        scale_x_continuous(breaks = c(0, 3, 8, 16)),
                        labs(x = "", y = "", fill = "Enzyme"), theme_bw(), theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                 panel.background = element_blank(), axis.text = element_text(color="black"), axis.line = element_line(colour = "black")))

# Functions to pivot dataframe or dataframe pair
pivot_taxonomy = function(df, cutoff=0) {
  df = df %>% 
    pivot_longer(
      cols = sample_list, 
      names_to = "sample",
      values_to = "tpm"
    ) %>% aggregate(tpm ~ sample + taxonomy, FUN = sum)
  df$sample = factor(df$sample, levels = sample_list)
  df = df[df$tpm > cutoff,]
  
  return(df)
}

pivot_annotation = function(df, cutoff=0) {
  df = df %>% 
    pivot_longer(
      cols = sample_list, 
      names_to = "sample",
      values_to = "tpm"
    ) %>% aggregate(tpm ~ sample + parent_bin + annotation, FUN = sum)
  df$sample = factor(df$sample, levels = sample_list)
  df = merge(x = df, y = parent_bin, by = "parent_bin", all.x = TRUE)
  df = df[df$tpm > cutoff,]
  
  return(df)
}

pivot_parent_bin = function(df_general, df_specific, cutoff=0) {
  df_general = df_general %>% 
    pivot_longer(
      cols = sample_list, 
      names_to = "sample",
      values_to = "tpm_total"
    ) %>% aggregate(tpm_total ~ sample + parent_bin, FUN = sum)
  df_general$sample = factor(df_general$sample, levels = sample_list)
  
  df_specific = df_specific %>% 
    pivot_longer(
      cols = sample_list, 
      names_to = "sample",
      values_to = "tpm_relevant"
    ) %>% aggregate(tpm_relevant ~ sample + parent_bin, FUN = sum)
  df_specific$sample = factor(df_specific$sample, levels = sample_list)
  
  df_general = merge(x = df_general, y = parent_bin, by = "parent_bin", all.x = TRUE)
  df_specific = merge(x = df_specific, y = df_general, by = c("sample","parent_bin"), all.x = TRUE)
  
  df_specific = df_specific[df_specific$tpm_total > cutoff,]
  df_specific$pct_relevant = (df_specific$tpm_relevant + 1) / (df_specific$tpm_total + 1)
  return(df_specific)
}

substrate_area_preprocess = function(df, filt, samp, substrate) {
  df_pivot = pivot_annotation(df[Reduce("|", lapply(filt, grepl, df$annotation)), ])
  df_pivot$annotation = ifelse(sub("\\_.*", "", df_pivot$annotation) != "S1", sub("\\_.*", "", df_pivot$annotation), df_pivot$annotation)
  df_pivot$annotation = sub("^.+\\+", "", df_pivot$annotation)
  df_pivot = aggregate(tpm ~ sample + annotation, df_pivot, sum, all.x=TRUE, drop = F)
  df_pivot = df_pivot[df_pivot$sample %in% samp,]
  df_pivot$sample <- as.numeric(as.character(mapvalues(df_pivot$sample, samp, c("0", "3", "8", "16"))))
  df_pivot$substrate = substrate
  return(df_pivot)
}

substrate_tax_area_preprocess = function(df, filt, samp, substrate) {
  df_pivot = pivot_taxonomy(df[Reduce("|", lapply(filt, grepl, df$annotation)), ])
  df_pivot = aggregate(tpm ~ sample + taxonomy, df_pivot, sum, all.x=TRUE, drop = F)
  df_pivot = df_pivot[df_pivot$sample %in% samp,]
  df_pivot$sample <- as.numeric(as.character(mapvalues(df_pivot$sample, samp, c("0", "3", "8", "16"))))
  df_pivot$substrate = substrate
  return(df_pivot)
}

annot_aggregate_area = function(df, samp, type="enzyme") {
  combined = list()
  if (type == "enzyme") {
    combined = do.call(rbind, list(
      substrate_area_preprocess(df, alginate_list, samp, "Alginate"),
      substrate_area_preprocess(df, agar_list, samp, "Agar"),
      substrate_area_preprocess(df, carrageenan_list, samp, "Carrageenan"),
      substrate_area_preprocess(df, fucoidan_list, samp, "Fucoidan")))
  } else if (type == "detox") {
    combined = do.call(rbind, list(
      substrate_area_preprocess(df, c("PGR"), samp, "Phloroglucinol"),
      substrate_area_preprocess(df, c("CarA", "CarB", "CarC", "CarD", "CarE"), samp, "Caffeate"),
      substrate_area_preprocess(df, c("Catalase"), samp, "ROS"),
      substrate_area_preprocess(df, c("Dehalogenase"), samp, "Halogenateds")))
  } else {
    combined = do.call(rbind, list(
      substrate_tax_area_preprocess(df, alginate_list, samp, "Alginate"),
      substrate_tax_area_preprocess(df, carrageenan_list, samp, "Carrageenan"),
      substrate_tax_area_preprocess(df, fucoidan_list, samp, "Fucoidan"),
      substrate_tax_area_preprocess(df, ulvan_list, samp, "Ulvan")))
  }
  
  return(aggregate(tpm ~ sample + substrate, combined, sum, all.x=TRUE, drop = F))
}

enzyme_area_figure = function(df, genes, samp) {
  df = substrate_area_preprocess(df, genes, samp, "")
  df[is.na(df)] = 0
  df = df %>% mutate(annotation = forcats::fct_reorder(.f = annotation, 
                                                       .x = tpm,
                                                       .fun = function(.x) sum(.x),
                                                       .desc = FALSE))
  print(df)
  
  positioning<-df%>%
    filter(sample==3)%>%              # Keep only 2017 value
    arrange(desc(annotation))%>%                # Inverse factor order (first is at the bottom of plot)
  mutate(                              # Create new column ypos and
      ypos=cumsum(tpm)    # fill with cumulative sum of invest for 2017
    )
  positioning$ypos = ifelse(positioning$tpm < -1, -1000, positioning$ypos - (positioning$tpm / 2))
  
  plot = ggplot(df, mapping = aes(x = sample, y = tpm, group = annotation, fill = annotation)) + area_plot_params + theme(legend.position="none") +
   geom_label(
      data=positioning, direction = "y",            # Different data than the main plot
      aes(y=ypos,label= annotation),   # y and lab in the aes() (change between labels)
      x=5, force = 0.5, min.segment.length = 0, size = 2                 # x is the same for all labels
    ) + theme(legend.position="none")
  return(plot)
}

tax_area_figure = function(df, genes, samp) {
  df = substrate_tax_area_preprocess(df, genes, samp, "")
  plot = ggplot(df, mapping = aes(x = sample, y = tpm, group = taxonomy, fill = taxonomy)) + area_plot_params + scale_fill_manual(values = taxonomy_coloring)
  return(plot)
}

# Order the x-axis samples and categories of protein regex
sample_list = c("FISH_1", "FISH_2", "FISH_3", "FISH_4", "FISH_5", "FISH_6", "T0_A", "T0_B", "T0_C", "T0_D", "M30_T03", "M30_T08", "M30_T16", "M36_T03", "M36_T08", "M36_T16", "M75_T03", "M75_T08", "M75_T16", "M78_T03", "M78_T08", "M78_T16")

alginate_list = c("PL6$", "PL6_", "PL7$", "PL7_", "PL17", "PL32", "PL34", "PL36", "PL38")
alginate_colors = c("PL6" = "#E5D8C7", "PL7" = "#9A7647", "PL17" = "#624C2D", "PL32" = "#7E613A", "PL34" = "#463620", "PL36" = "#B89465", "PL38"  = "#D2BB9D")

carrageenan_list = c("S1_7$", "S1_19", "GH16_17" , "GH127", "GH129", "GH82")
carrageenan_colors = c("S1_7" = "#9A133DFF", "S1_19" = "#B93961FF", "GH129" = "#D8527CFF", "GH16" = "#F9E0E8FF", "GH82" = "#F9B4C9FF", "GH127" = "#F28AAAFF" )

fucoidan_list = c("S1_15",	"S1_16",	"S1_17", "S1_25",
                  "GH29",	"GH95","GH168",	"GH141", "GH107")
fucoidan_colors = c("S1_15" = "#ECEFF1FF", "S1_16" = "#CFD8DCFF", "S1_17" = "#B0BEC5FF", "S1_25" = "#37474FFF",
                    "GH29" = "#90A4AEFF", "GH95" = "#78909CFF", "GH107" = "#607D8BFF", "GH168" = "#546E7AFF", "GH141" = "#455A64FF")

substrate_colors = c("Agar" = "red", "Alginate" = "#7E613A", "Carrageenan" = "lightpink", "Fucoidan" = "#CFD8DCFF", "Ulvan" = "green", "Phloroglucinol" = "grey40","Caffeate" = "black", "ROS" = "red", "Halogenateds"="blue")

ulvan_list = c("PL24", "PL28", "PL25", "PL40")
agar_list = c("GH117", "GH118", "GH50", "GH86", "GH96", "GH16_14", "GH16_15", "GH16_16")


# Import Expression data, protein and MAG descriptions
mapping_rate = read.csv("data/mapping_rate.tab", sep = "\t")
tax_table = read.csv("data/bin_taxonomy.tab", sep = "\t")
parent_bin = unique.data.frame(tax_table[c("taxonomy", "parent_bin")])

rsem = left_join(read.csv("data/binned_expression_fixed.tab", sep = "\t"),
                 read.csv("data/cazyme.tab", sep = "\t"), by= "protein") %>%
       left_join(., read.csv("data/sulfatase.tab", sep = "\t"), by = "protein") %>%
       left_join(., read.csv("data/camper.tab", sep = "\t"), by = "protein") %>%
       left_join(., read.csv("data/protein_bins.tab", sep = "\t"), by = "protein") %>%
       left_join(., read.csv("data/COGS.tab", sep = "\t"), by = "protein") %>%
       left_join(., read.csv("data/cog_cats.tab", sep = "\t"), by = "COG") %>%
       left_join(., read.csv("data/detox.tab", sep = "\t"), by = "protein") %>%
       left_join(., tax_table, by = "bin")
  

#ECEFF1FF, #CFD8DCFF, #B0BEC5FF, #90A4AEFF, #78909CFF, #607D8BFF, #546E7AFF, #455A64FF, #37474FFF, #263238FF
# Normalize TPMs for mapping rates (only 20-60% of reads are MAG-resolved)
rsem_norm = rsem
for (samp in sample_list) {
  rsem_norm[[samp]] = rsem_norm[[samp]] * (dplyr::filter(mapping_rate, sample==samp)$map_rate * 0.01)
  rsem_norm[[samp]] = (rsem_norm[[samp]] + 0.00000001)
}

rsem_norm$T0_AVG = rowMeans(rsem_norm[,c("T0_A","T0_B","T0_C","T0_D")])
rsem_norm = rsem_norm[, !(names(rsem_norm) %in% c("T0_A","T0_B","T0_C","T0_D"))]
sample_list = c("FISH_1", "FISH_2", "FISH_3", "FISH_4", "FISH_5", "FISH_6", "T0_AVG", "M30_T03", "M30_T08", "M30_T16", "M36_T03", "M36_T08", "M36_T16", "M75_T03", "M75_T08", "M75_T16", "M78_T03", "M78_T08", "M78_T16")


# Drop poorly classified MAGs and standardize annotation column
rsem_norm = rsem_norm[rsem_norm$taxonomy != "Unclassified",]
rsem_norm$annotation = ifelse(!is.na(rsem_norm$sulfatase), rsem_norm$sulfatase, ifelse(!is.na(rsem_norm$camper), rsem_norm$camper, ifelse(!is.na(rsem_norm$CAZyme), rsem_norm$CAZyme, rsem_norm$detox)))

###### AREA PLOTS

### DETOX

sample_m30 = c("T0_AVG", "M30_T03", "M30_T08", "M30_T16")
sample_m36 = c("T0_AVG", "M36_T03", "M36_T08", "M36_T16")
sample_m75 = c("T0_AVG", "M75_T03", "M75_T08", "M75_T16")
sample_m78 = c("T0_AVG", "M78_T03", "M78_T08", "M78_T16")

m30 = ggplot(annot_aggregate_area(rsem_norm, sample_m30, type = "detox"), mapping = aes(x = sample, y = tpm, fill = substrate)) +
  geom_area(stat="identity", position = "stack") + area_plot_params + labs(fill = "Substrate") + #ylim(0, 2500) + 
  scale_fill_manual(values = substrate_colors)

m36 = ggplot(annot_aggregate_area(rsem_norm, sample_m36, type = "detox"), mapping = aes(x = sample, y = tpm, fill = substrate)) +
  geom_area(stat="identity", position = "stack") + area_plot_params + labs(fill = "Substrate") + #ylim(0, 2500) + 
  scale_fill_manual(values = substrate_colors) + theme(legend.position="none")

m75 = ggplot(annot_aggregate_area(rsem_norm, sample_m75, type = "detox"), mapping = aes(x = sample, y = tpm, fill = substrate)) +
  geom_area(stat="identity", position = "stack") + area_plot_params + labs(fill = "Substrate") + #ylim(0, 2500) + 
  scale_fill_manual(values = substrate_colors) + theme(legend.position="none")

m78 = ggplot(annot_aggregate_area(rsem_norm, sample_m78, type = "detox"), mapping = aes(x = sample, y = tpm, fill = substrate)) +
  geom_area(stat="identity", position = "stack") + area_plot_params + labs(fill = "Substrate") + #ylim(0, 2500) + 
  scale_fill_manual(values = substrate_colors) + theme(legend.position="none")

gridded = plot_grid(
  plot_grid(
    m30 + theme(legend.position="none") +  ylab("TPM"), m36, m75, m78,
    nrow = 1), 
  
  get_legend(m30 + theme(legend.box.margin = margin(0, 0, 0, 12))), rel_widths = c(3, .6))

gridded

caff30_tax = tax_area_figure(rsem_norm, c("CarA", "CarB", "CarC", "CarD", "CarE"), sample_m30) +ylim(0, NA) 
caff36_tax  = tax_area_figure(rsem_norm, c("CarA", "CarB", "CarC", "CarD", "CarE"), sample_m36) +ylim(0, NA)
caff75_tax  = tax_area_figure(rsem_norm, c("CarA", "CarB", "CarC", "CarD", "CarE"), sample_m75) +ylim(0, NA) 
caff78_tax  = tax_area_figure(rsem_norm, c("CarA", "CarB", "CarC", "CarD", "CarE"), sample_m78) + ylim(0, NA)

pgr30_tax  = tax_area_figure(rsem_norm, c("PGR"), sample_m30) + ylim(0, NA) 
pgr36_tax  = tax_area_figure(rsem_norm, c("PGR"), sample_m36) + ylim(0, NA) 
pgr75_tax  = tax_area_figure(rsem_norm, c("PGR"), sample_m75) + ylim(0, NA)  
pgr78_tax  = tax_area_figure(rsem_norm, c("PGR"), sample_m78) + ylim(0, NA)

ros30_tax  = tax_area_figure(rsem_norm, c("Catalase"), sample_m30) + ylim(0, NA)  
ros36_tax  = tax_area_figure(rsem_norm, c("Catalase"), sample_m36) + ylim(0, NA) 
ros75_tax  = tax_area_figure(rsem_norm, c("Catalase"), sample_m75) + ylim(0, NA)  
ros78_tax  = tax_area_figure(rsem_norm, c("Catalase"), sample_m78) + ylim(0, NA)  

big_gridded = plot_grid(plot_grid(
  m30 + theme(legend.position="none") + ylab("TPM"), m36, m75, m78,
  caff30_tax + ylab("TPM")+ theme(legend.position="none"), caff36_tax+ theme(legend.position="none"), caff75_tax+ theme(legend.position="none"), caff78_tax+ theme(legend.position="none"),
  pgr30_tax + ylab("TPM")+ theme(legend.position="none"), pgr36_tax+ theme(legend.position="none"), pgr75_tax+ theme(legend.position="none"), pgr78_tax+ theme(legend.position="none"),
  ros30_tax + ylab("TPM")+ theme(legend.position="none"), ros36_tax+ theme(legend.position="none"), ros75_tax+ theme(legend.position="none"), ros78_tax+ theme(legend.position="none"), nrow = 4, labels = 
    c("A", "", "", "", "B","", "", "",  "C", "", "", "", "D","", "", "")), plot_grid(get_legend(m30 + theme(legend.box.margin = margin(0, 0, 0, 12))), get_legend(caff30_tax + theme(legend.box.margin = margin(0, 0, 0, 12))), nrow = 2, rel_heights = c(0.25, 0.75)),
  rel_widths = c(3, .6), nrow = 1)

svg("combined_detox_tax_distributions.svg", width = 10, height = 6)
big_gridded
dev.off()

#### CARBO

sample_m30 = c("T0_AVG", "M30_T03", "M30_T08", "M30_T16")
sample_m36 = c("T0_AVG", "M36_T03", "M36_T08", "M36_T16")
sample_m75 = c("T0_AVG", "M75_T03", "M75_T08", "M75_T16")
sample_m78 = c("T0_AVG", "M78_T03", "M78_T08", "M78_T16")

m30 = ggplot(annot_aggregate_area(rsem_norm, sample_m30), mapping = aes(x = sample, y = tpm, fill = substrate)) +
  geom_area(stat="identity", position = "stack") + area_plot_params + labs(fill = "Substrate") + #ylim(0, 2500) + 
  scale_fill_manual(values = substrate_colors)

m36 = ggplot(annot_aggregate_area(rsem_norm, sample_m36), mapping = aes(x = sample, y = tpm, fill = substrate)) +
  geom_area(stat="identity", position = "stack") + area_plot_params + labs(fill = "Substrate") + #ylim(0, 2500) + 
  scale_fill_manual(values = substrate_colors) + theme(legend.position="none")
  
m75 = ggplot(annot_aggregate_area(rsem_norm, sample_m75), mapping = aes(x = sample, y = tpm, fill = substrate)) +
  geom_area(stat="identity", position = "stack") + area_plot_params + labs(fill = "Substrate") + #ylim(0, 2500) + 
  scale_fill_manual(values = substrate_colors) + theme(legend.position="none")
  
m78 = ggplot(annot_aggregate_area(rsem_norm, sample_m78), mapping = aes(x = sample, y = tpm, fill = substrate)) +
  geom_area(stat="identity", position = "stack") + area_plot_params + labs(fill = "Substrate") + #ylim(0, 2500) + 
  scale_fill_manual(values = substrate_colors) + theme(legend.position="none")

gridded = plot_grid(
  plot_grid(
    m30 + theme(legend.position="none") +  ylab("TPM"), m36, m75, m78,
    nrow = 1), 
  
  get_legend(m30 + theme(legend.box.margin = margin(0, 0, 0, 12))), rel_widths = c(3, .6))

gridded

alg30 = enzyme_area_figure(rsem_norm, alginate_list, sample_m30) + scale_fill_manual(values = alginate_colors)
alg36 = enzyme_area_figure(rsem_norm, alginate_list, sample_m36) + scale_fill_manual(values = alginate_colors)
alg75 = enzyme_area_figure(rsem_norm, alginate_list, sample_m75) + scale_fill_manual(values = alginate_colors)
alg78 = enzyme_area_figure(rsem_norm, alginate_list, sample_m78) + scale_fill_manual(values = alginate_colors)

carr30 = enzyme_area_figure(rsem_norm, carrageenan_list, sample_m30) + scale_fill_manual(values = carrageenan_colors)
carr36 = enzyme_area_figure(rsem_norm, carrageenan_list, sample_m36) + scale_fill_manual(values = carrageenan_colors)
carr75 = enzyme_area_figure(rsem_norm, carrageenan_list, sample_m75) + scale_fill_manual(values = carrageenan_colors)
carr78 = enzyme_area_figure(rsem_norm, carrageenan_list, sample_m78) + scale_fill_manual(values = carrageenan_colors)

fuc30 = enzyme_area_figure(rsem_norm, fucoidan_list, sample_m30) + scale_fill_manual(values = fucoidan_colors)
fuc36 = enzyme_area_figure(rsem_norm, fucoidan_list, sample_m36) + scale_fill_manual(values = fucoidan_colors)
fuc75 = enzyme_area_figure(rsem_norm, fucoidan_list, sample_m75) + scale_fill_manual(values = fucoidan_colors)
fuc78 = enzyme_area_figure(rsem_norm, fucoidan_list, sample_m78) + scale_fill_manual(values = fucoidan_colors)


big_gridded = plot_grid(plot_grid(
  m30 + theme(legend.position="none") + ylab("TPM"), m36, m75, m78,
  alg30 + ylab("TPM"), alg36, alg75, alg78,
  carr30 + ylab("TPM"), carr36, carr75, carr78,
  fuc30 + ylab("TPM"), fuc36, fuc75, fuc78, nrow = 4, labels = 
    c("A", "", "", "", "B","", "", "",  "C", "", "", "", "D","", "", "")), plot_grid(get_legend(m30 +  theme(legend.justification = "top")), theme_nothing(), theme_nothing(), theme_nothing(), nrow = 4, rel_heights = c(0.25, 0.25, 0.25, 0.25)),
rel_widths = c(3, .6), nrow = 1)
big_gridded

svg("combined_wagar_label_enzyme_distributions.svg", width = 10, height = 6)
big_gridded
dev.off()

alg30 = enzyme_area_figure(rsem_norm, alginate_list, sample_m30) +ylim(0, NA) +  scale_fill_manual(values = alginate_colors)
alg36 = enzyme_area_figure(rsem_norm, alginate_list, sample_m36) +ylim(0, NA) +  scale_fill_manual(values = alginate_colors)
alg75 = enzyme_area_figure(rsem_norm, alginate_list, sample_m75) +ylim(0, NA) +  scale_fill_manual(values = alginate_colors)
alg78 = enzyme_area_figure(rsem_norm, alginate_list, sample_m78) + ylim(0, NA) + scale_fill_manual(values = alginate_colors)

carr30 = enzyme_area_figure(rsem_norm, carrageenan_list, sample_m30) + ylim(0, NA) + scale_fill_manual(values = carrageenan_colors)
carr36 = enzyme_area_figure(rsem_norm, carrageenan_list, sample_m36) + ylim(0, NA) + scale_fill_manual(values = carrageenan_colors)
carr75 = enzyme_area_figure(rsem_norm, carrageenan_list, sample_m75) + ylim(0, NA) +  scale_fill_manual(values = carrageenan_colors)
carr78 = enzyme_area_figure(rsem_norm, carrageenan_list, sample_m78) + ylim(0, NA) +  scale_fill_manual(values = carrageenan_colors)

fuc30 = enzyme_area_figure(rsem_norm, fucoidan_list, sample_m30) + ylim(0, NA) +  scale_fill_manual(values = fucoidan_colors)
fuc36 = enzyme_area_figure(rsem_norm, fucoidan_list, sample_m36) + ylim(0, NA) + scale_fill_manual(values = fucoidan_colors)
fuc75 = enzyme_area_figure(rsem_norm, fucoidan_list, sample_m75) + ylim(0, NA) +  scale_fill_manual(values = fucoidan_colors)
fuc78 = enzyme_area_figure(rsem_norm, fucoidan_list, sample_m78) + ylim(0, NA) +  scale_fill_manual(values = fucoidan_colors)


big_gridded = plot_grid(plot_grid(
  m30 + theme(legend.position="none") + ylab("TPM"), m36, m75, m78,
  alg30 + ylab("TPM"), alg36, alg75, alg78,
  carr30 + ylab("TPM"), carr36, carr75, carr78,
  fuc30 + ylab("TPM"), fuc36, fuc75, fuc78, nrow = 4, labels = 
    c("A", "", "", "", "B","", "", "",  "C", "", "", "", "D","", "", "")), plot_grid(get_legend(m30 + theme(legend.box.margin = margin(0, 0, 0, 12))), get_legend(alg30 + theme(legend.box.margin = margin(0, 0, 0, 12))), theme_nothing(), theme_nothing(), nrow = 4),
  rel_widths = c(3, .6), nrow = 1)

big_gridded

svg("combined_enzyme_distributions_nolim.svg", width = 10, height = 6)
big_gridded
dev.off()

alg30_tax = tax_area_figure(rsem_norm, alginate_list, sample_m30) +ylim(0, NA) 
alg36_tax  = tax_area_figure(rsem_norm, alginate_list, sample_m36) +ylim(0, NA)
alg75_tax  = tax_area_figure(rsem_norm, alginate_list, sample_m75) +ylim(0, NA) 
alg78_tax  = tax_area_figure(rsem_norm, alginate_list, sample_m78) + ylim(0, NA)

carr30_tax  = tax_area_figure(rsem_norm, carrageenan_list, sample_m30) + ylim(0, NA) 
carr36_tax  = tax_area_figure(rsem_norm, carrageenan_list, sample_m36) + ylim(0, NA) 
carr75_tax  = tax_area_figure(rsem_norm, carrageenan_list, sample_m75) + ylim(0, NA)  
carr78_tax  = tax_area_figure(rsem_norm, carrageenan_list, sample_m78) + ylim(0, NA)
s
fuc30_tax  = tax_area_figure(rsem_norm, fucoidan_list, sample_m30) + ylim(0, NA)  
fuc36_tax  = tax_area_figure(rsem_norm, fucoidan_list, sample_m36) + ylim(0, NA) 
fuc75_tax  = tax_area_figure(rsem_norm, fucoidan_list, sample_m75) + ylim(0, NA)  
fuc78_tax  = tax_area_figure(rsem_norm, fucoidan_list, sample_m78) + ylim(0, NA)  

big_gridded = plot_grid(plot_grid(
  m30 + theme(legend.position="none") + ylab("TPM"), m36, m75, m78,
  alg30_tax + ylab("TPM")+ theme(legend.position="none"), alg36_tax+ theme(legend.position="none"), alg75_tax+ theme(legend.position="none"), alg78_tax+ theme(legend.position="none"),
  carr30_tax + ylab("TPM")+ theme(legend.position="none"), carr36_tax+ theme(legend.position="none"), carr75_tax+ theme(legend.position="none"), carr78_tax+ theme(legend.position="none"),
  fuc30_tax + ylab("TPM")+ theme(legend.position="none"), fuc36_tax+ theme(legend.position="none"), fuc75_tax+ theme(legend.position="none"), fuc78_tax+ theme(legend.position="none"), nrow = 4, labels = 
    c("A", "", "", "", "B","", "", "",  "C", "", "", "", "D","", "", "")), plot_grid(get_legend(m30 + theme(legend.box.margin = margin(0, 0, 0, 12))), get_legend(alg30_tax + theme(legend.box.margin = margin(0, 0, 0, 12))), nrow = 2, rel_heights = c(0.25, 0.75)),
  rel_widths = c(3, .6), nrow = 1)

big_gridded
tax_area_figure(rsem_norm, fucoidan_list, sample_m78) + ylim(0, NA)

svg("combined_enzyme_tax_distributions_nolim.svg", width = 10, height = 6)
big_gridded
dev.off()






##### Expression based on metaSPAdes assembled metatranscripts
assembled_rsem = read.csv("data/assembled_expression.tab", sep = "\t")  %>% 
  pivot_longer(
    cols = c("CAZyme", "Sulfatase"),
    names_to = "Annotation",
    values_to = "Expression"
  )
assembled_rsem$sample = factor(assembled_rsem$sample, levels = sample_list)

svg("figures/metaT_assembled_expression.svg", width=10, height=6)
ggplot(assembled_rsem, aes(x = sample, y = Expression, fill=Annotation, group=Annotation)) +
  geom_bar(stat="identity", width=.5, position = "dodge") + theme_bw() +
  scale_fill_manual(name = "Gene Annotation", values = c("CAZyme" = "black", "Sulfatase" = "gray")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")) + 
  ylab("metaSPAdes Assembled Sample TPM") + xlab("Metatranscriptomic Sample")
dev.off()

scale_fill_manual(name = "MAG Taxonomy", values = taxonomy_coloring) +
  
  geom_text(data = test %>% filter(sample == 0),
            aes(x = sample - 0.5, y = tpm, 
                label = annotation, color = annotation),
            position = position_stack(vjust = 0.5))


####

rsem_norm_parent = rsem_norm[!is.na(rsem_norm$CAZyme) | !is.na(rsem_norm$sulfatase),]
rsem_norm_parent = rsem_norm_parent %>%  group_by(parent_bin) %>%
  dplyr::summarise(across(where(is.numeric), sum), .groups = 'drop')


rsem_norm_COG = rsem_norm[!is.na(rsem_norm$cat),]
rsem_norm_COG = rsem_norm_COG %>%  group_by(cat) %>%
  dplyr::summarise(across(where(is.numeric), sum), .groups = 'drop')

rsem_norm_parent = aggregate(.~parent_bin, rsem_norm_parent, sum)

summarise(across(where(is.numeric), sum), .groups = 'drop')

## Total Expression Figure
total_expression_plot = ggplot(pivot_taxomony(rsem_norm), 
                               aes(fill=taxonomy, y=tpm, x=sample)) +
  bar_graph_params + ylab("Genome-resolved TPM")

svg("figures/metaT_resolved_total.svg", width=10, height=6)
total_expression_plot
dev.off()

## Bulk CAZyme Expression
{
  total_cazy_expression_plot = ggplot(pivot_taxomony(rsem_norm[is.na(rsem_norm$sulfatase) & !is.na(rsem_norm$CAZyme),]),
                                      aes(fill=taxonomy, y=tpm, x=sample)) + 
    bar_graph_params + ylab("Genome-resolved CAZyme TPM")
  
  svg("figures/metaT_resolved_total_CAZyme.svg", width=10, height=6)
  total_cazy_expression_plot 
  dev.off()
  
  cazy_expr = pivot_parent_bin(rsem_norm, rsem_norm[is.na(rsem_norm$sulfatase) & !is.na(rsem_norm$CAZyme), ], 100)
  
  pct_cazy_expression_plot = ggplot(pivot_parent_bin(rsem_norm, rsem_norm[is.na(rsem_norm$sulfatase) & !is.na(rsem_norm$CAZyme), ], 5000),
                                    aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
    beeswarm_params + ylab("CAZymes as % of MAG-specific Expression")  + 
    scale_y_continuous(labels = scales::percent)
  
  svg("figures/metaT_resolved_CAZyme_pct.svg", width=10, height=6)
  pct_cazy_expression_plot
  dev.off()
}

## Sulfatases
{
  total_sulf_expression_plot = ggplot(pivot_taxomony(rsem_norm[!is.na(rsem_norm$sulfatase),]),
                                      aes(fill=taxonomy, y=tpm, x=sample)) + 
    bar_graph_params + ylab("Genome-resolved Sulfatase TPM")
  
  svg("figures/metaT_resolved_total_sulf.svg", width=10, height=6)
  total_sulf_expression_plot
  dev.off()
  
  sulf_expr = pivot_parent_bin(rsem_norm, rsem_norm[!is.na(rsem_norm$sulfatase), ], 100)
  detox_expr = pivot_parent_bin(rsem_norm, rsem_norm[!is.na(rsem_norm$camper) | !is.na(rsem_norm$detox), ], 100)
  
  
  pct_sulfatase_expression_plot = ggplot(pivot_parent_bin(rsem_norm, rsem_norm[!is.na(rsem_norm$sulfatase), ], 5000),
                                         aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
    beeswarm_params + ylab("Sulfatases as % of MAG-specific Expression")  + 
    scale_y_continuous(labels = scales::percent)
  
  svg("figures/metaT_resolved_sulf_pct.svg", width=10, height=6)
  pct_sulfatase_expression_plot
  dev.off()
}

## Polyphenol Degradation
{
  total_camper_expression_plot = ggplot(pivot_taxomony(rsem_norm[!is.na(rsem_norm$sulfatase),]),
                                        aes(fill=taxonomy, y=tpm, x=sample)) + 
    bar_graph_params + ylab("Genome-resolved Polyphenol-degrading TPM")
  
  svg("figures/metaT_resolved_total_camper.svg", width=10, height=6)
  total_camper_expression_plot
  dev.off()
}

## Top Individual Gene Expression - Alginate
{
  top_alginate_expression_plot = ggplot(pivot_annotation(rsem_norm[Reduce("|", lapply(alginate_list, grepl, rsem_norm$annotation)), ], 10),
                                        aes(x = sample, y = tpm, color=taxonomy, group=parent_bin)) +
    beeswarm_params + ylab("Alginate Enzyme TPM")
  
  svg("figures/metaT_resolved_total_alginate.svg", width=10, height=6)
  top_alginate_expression_plot 
  dev.off()
  
  pct_alginate_expression_plot = ggplot(pivot_parent_bin(rsem_norm, rsem_norm[Reduce("|", lapply(alginate_list, grepl, rsem_norm$annotation)), ], 5000),
                                        aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
    beeswarm_params + ylab("Alginate-degradation as % of MAG-specific Expression")  + 
    scale_y_continuous(labels = scales::percent)
  
  svg("figures/metaT_resolved_alginate_pct.svg", width=10, height=6)
  pct_alginate_expression_plot
  dev.off()
}

## Top Individual Gene Expression - Carrageenan
{
  top_carrageenan_expression_plot = ggplot(pivot_annotation(rsem_norm[Reduce("|", lapply(carrageenan_list, grepl, rsem_norm$annotation)), ], 10),
                                           aes(x = sample, y = tpm, color=taxonomy, group=parent_bin)) +
    beeswarm_params + ylab("Carrageenan Enzyme TPM")
  
  svg("figures/metaT_resolved_total_carra.svg", width=10, height=6)
  top_carrageenan_expression_plot 
  dev.off()
  
  pct_carrageenan_expression_plot = ggplot(pivot_parent_bin(rsem_norm, rsem_norm[Reduce("|", lapply(carrageenan_list, grepl, rsem_norm$annotation)), ], 5000),
                                           aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
    beeswarm_params + ylab("Carrageenan-degradation as % of MAG-specific Expression")  + 
    scale_y_continuous(labels = scales::percent)
  
  svg("figures/metaT_resolved_carra_pct.svg", width=10, height=6)
  pct_carrageenan_expression_plot
  dev.off()
}

##  Top Individual Gene Expression - Fucoidan
{
  top_fucoidan_expression_plot = ggplot(pivot_annotation(rsem_norm[Reduce("|", lapply(fucoidan_list, grepl, rsem_norm$annotation)), ], 10),
                                        aes(x = sample, y = tpm, color=taxonomy, group=parent_bin)) +
    beeswarm_params + ylab("Fucoidan Enzyme TPM")
  
  svg("figures/metaT_resolved_total_fcsp.svg", width=10, height=6)
  top_fucoidan_expression_plot 
  dev.off()
  
  pct_fucoidan_expression_plot = ggplot(pivot_parent_bin(rsem_norm, rsem_norm[Reduce("|", lapply(fucoidan_list, grepl, rsem_norm$annotation)), ], 5000),
                                        aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
    beeswarm_params + ylab("Fucoidan-degradation as % of MAG-specific Expression") + 
    scale_y_continuous(labels = scales::percent)
  
  svg("figures/metaT_resolved_fcsp_pct.svg", width=10, height=6)
  pct_fucoidan_expression_plot
  dev.off()
}



##########

cog_plot_params = list(geom_area(stat="identity", position = "fill", color = "black"),
                        scale_x_continuous(breaks = c(0, 3, 8, 16)), scale_y_continuous(labels = scales::percent_format()),
                        labs(x = "", y = "", fill = "COG Category"), theme_bw(), theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                 panel.background = element_blank(), axis.text = element_text(color="black"),  axis.line = element_line(colour = "black")))

cogs = read.csv("cogs.tab", sep = "\t")
cog_level_list = read.csv("cogs_levels.tab", sep = "\t", header = F)$V1
cogs$name = factor(cogs$name, levels = cog_level_list)
cogs = cogs[!(cogs$cat %in% c("A", "W", "X", "Z", "S", "Q")),]


cogs_m30 = cogs[cogs$microcosm %in% c("m30", "T0"),]
cogs_m30$timepoint = as.numeric(cogs_m30$timepoint)
m30_plot = ggplot(cogs_m30, mapping = aes(x = timepoint, y = tpm, fill = name)) + cog_plot_params

cogs_m36 = cogs[cogs$microcosm %in% c("m36", "T0"),]
cogs_m36$timepoint = as.numeric(cogs_m36$timepoint)
m36_plot = ggplot(cogs_m36, mapping = aes(x = timepoint, y = tpm, fill = name)) + cog_plot_params

cogs_m75 = cogs[cogs$microcosm %in% c("m75", "T0"),]
cogs_m75$timepoint = as.numeric(cogs_m75$timepoint)
m75_plot = ggplot(cogs_m75, mapping = aes(x = timepoint, y = tpm, fill = name)) + cog_plot_params

cogs_m78 = cogs[cogs$microcosm %in% c("m78", "T0"),]
cogs_m78$timepoint = as.numeric(cogs_m78$timepoint)
m78_plot = ggplot(cogs_m78, mapping = aes(x = timepoint, y = tpm, fill = name)) + cog_plot_params

cog_plot = plot_grid(m30_plot + theme(legend.position="none"), m36_plot + theme(legend.position="none"), m75_plot + theme(legend.position="none"), m78_plot + theme(legend.position="none"), get_legend(m30_plot), nrow = 1, rel_widths = c(0.5, 0.5, 0.5, 0.5, 1))

svg("cog_plot.svg", width=12, height=6)
cog_plot
dev.off()



max_pct = read.table("max_pct_1k.tab", sep = "\t", header = T)

caz = ggplot(max_pct[!is.na(max_pct$max_cazy_pct),], aes(x=pct_cazy_genome, y=max_cazy_pct)) + geom_point(size = 3, aes(color=taxonomic.group)) + scale_color_manual(values = taxonomy_coloring) + theme_bw()  + 
                            scale_y_continuous(labels = scales::percent)  + 
                            scale_x_continuous(labels = scales::percent) + guides(color = guide_legend(override.aes = list(size = 5)))  +
  guides(color=guide_legend(title="MAG Taxonomy", size = 14)) + theme(
    axis.text = element_text(color="black", size = 14),
    axis.ticks = element_line(color = "black"),
    axis.title=element_text(size=14)
                            ) + xlab("CAZymes as % of Genome") + ylab("Maximum CAZymes as % of Expression")
                          

sulf = ggplot(max_pct[!is.na(max_pct$max_sulf_pct),], aes(x=pct_sulf_genome, y=max_sulf_pct)) + geom_point(size = 3, aes(color=taxonomic.group)) + scale_color_manual(values = taxonomy_coloring) + theme_bw()  + 
                             scale_y_continuous(labels = scales::percent)  + 
                             scale_x_continuous(labels = scales::percent) + guides(color = guide_legend(override.aes = list(size = 5)))  +
  guides(color=guide_legend(title="MAG Taxonomy", size = 14)) + theme(
    axis.text = element_text(color="black", size = 14),
    axis.ticks = element_line(color = "black"),
    axis.title=element_text(size=14)
                             ) + xlab("Sulfatases as % of Genome") + ylab("Maximum Sulfatases as % of Expression")
                           

detox = ggplot(max_pct[!is.na(max_pct$max_detox_pct),], aes(x=pct_detox_genome, y=max_detox_pct)) + geom_point(size = 3, aes(color=taxonomic.group)) + scale_color_manual(values = taxonomy_coloring) + theme_bw()  + 
  scale_y_continuous(labels = scales::percent)  + 
  scale_x_continuous(labels = scales::percent) + guides(color = guide_legend(override.aes = list(size = 7))) +
  guides(color=guide_legend(title="MAG Taxonomy", size = 14)) + theme(
    axis.text = element_text(color="black", size = 14),
    axis.ticks = element_line(color = "black"),
    axis.title=element_text(size=14)
  ) + xlab("Detoxification Genes as % of Genome") + ylab("Maximum Detoxification Genes as % of Expression")

svg("caz_plot.svg", width=12, height=6)
caz
dev.off()

svg("sulf_plot.svg", width=12, height=6)
sulf
dev.off()

svg("detox_plot.svg", width=12, height=6)
detox
dev.off()

# stat_poly_line(se = F, color = "black", linewidth = 0.5) +
# stat_poly_eq(use_label(c("eq", "R2"))) + 

# Change the point size, and shape
ggplot(mtcars, aes(x=wt, y=mpg)) +
  geom_point(size=2, shape=23)