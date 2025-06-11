# Get visualization packages
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('ggbeeswarm')) install.packages('ggbeeswarm'); library('ggbeeswarm')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('ggtext')) install.packages('ggtext'); library('ggtext')
if (!require('ggnewscale')) install.packages('ggnewscale'); library('ggnewscale')


# Get data manipulation packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')

df = read.csv("C:/Users/aaron/Documents/krumbs_taxonomy.txt",sep="\t",stringsAsFactors = F)
df = df %>% pivot_longer(cols=c("Gammaproteobacteria", "Desulfobacterota", "Bacillota", "Alphaproteobacteria", "Bacteroidota", "Spirochaetota", "Verrucomicrobiota", "Fusobacterota", "Lentisphaerota", "Tenericutes"),
                    names_to='taxonomy',
                    values_to='abundance') 

taxa_order = c("Alphaproteobacteria", "Bacillota", "Bacteroidota", "Desulfobacterota", "Fusobacterota",
               "Gammaproteobacteria","Lentisphaerota", "Spirochaetota", "Tenericutes", "Verrucomicrobiota")
df$taxonomy = factor(df$taxonomy, levels = taxa_order)

#time_order = c("Fish", "T0", "T3", "T8", "T16")
#df$timepoint = factor(df$timepoint, levels = time_order)

#df_fishone = df[df$sample == "FISH_1_DNA",]

generate_cosm_area_plot = function(df_plot) {
  df_plot$timepoint = as.numeric(df_plot$timepoint)
  return(
         ggplot(df_plot, aes(x = timepoint, y = abundance, fill = taxonomy)) +
           geom_area() + 
           scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) + 
           scale_x_continuous(breaks = c(0, 3, 8, 16)) +
           labs(x = "", y = "", fill = "Read Taxonomy") + scale_fill_manual(labels = 
                                                                              c("Alphaproteobacteria", "Bacillota", "Bacteroidota", "Desulfovibrionales", "Fusobacterota",
                                                                                "Gammaproteobacteria","Lentisphaerota", "Spirochaetota", "Tenericutes", "Verrucomicrobiota"
                                                                              ), values = 
                                                                              c( color_alpha,  color_firm,  color_bacteroid,
                                                                                 color_desulfo,   
                                                                                 color_fuso, color_gamma, color_weird, color_spiro, color_tener,
                                                                                 color_pvc)) +
           theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.text = element_text(color="black"), axis.line = element_line(colour = "black")))
}

generate_fish_barplot = function(df_plot) {
    return(
      ggplot(df_plot, aes(x = sample, y = abundance, fill = taxonomy)) +
        geom_bar(stat = "identity") + 
        scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) + 
        labs(x = "", y = "", fill = "Read Taxonomy") + scale_fill_manual(labels = 
                                                                           c("Alphaproteobacteria", "Bacillota", "Bacteroidota", "Desulfovibrionales", "Fusobacterota",
                                                                             "Gammaproteobacteria","Lentisphaerota", "Spirochaetota", "Tenericutes", "Verrucomicrobiota"
                                                                           ), values = 
                                                                           c( color_alpha,  color_firm,  color_bacteroid,
                                                                              color_desulfo,   
                                                                              color_fuso, color_gamma, color_weird, color_spiro, color_tener,
                                                                              color_pvc)) +
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.text = element_text(color="black"), axis.line = element_line(colour = "black")))
  }

{
  color_alpha = "#006F6C"
  color_gamma = "#8453b5"
  color_bacteroid = "#006DDB"
  color_desulfo = "#22CF22"
  color_pvc = "#d5b2f7"
  color_spiro = "#ffdf4d"
  color_fuso = "#A7EAF7"
  color_firm = "#BC0000"
  color_tener = "#252525"
  color_weird = "#db6d00"
  color_euk = "#696969"
  color_other = "#D3D3D3"
}

one_mG = generate_cosm_area_plot(df[df$category %in% c("m30", "T0_AVG") & df$sequencing == "Metagenomic",])
one_mT = generate_cosm_area_plot(df[df$category %in% c("m30", "T0_AVG") & df$sequencing == "Metatranscriptomic",])
one_s = generate_cosm_area_plot(df[df$category %in% c("m30", "T0_AVG") & df$sequencing == "16S",])
two_mG = generate_cosm_area_plot(df[df$category %in% c("m36", "T0_AVG") & df$sequencing == "Metagenomic",])
two_mT = generate_cosm_area_plot(df[df$category %in% c("m36", "T0_AVG") & df$sequencing == "Metatranscriptomic",])
two_s = generate_cosm_area_plot(df[df$category %in% c("m36", "T0_AVG") & df$sequencing == "16S",])
three_mG = generate_cosm_area_plot(df[df$category %in% c("m75", "T0_AVG") & df$sequencing == "Metagenomic",])
three_mT = generate_cosm_area_plot(df[df$category %in% c("m75", "T0_AVG") & df$sequencing == "Metatranscriptomic",])
three_s = generate_cosm_area_plot(df[df$category %in% c("m75", "T0_AVG") & df$sequencing == "16S",])
four_mG = generate_cosm_area_plot(df[df$category %in% c("m78", "T0_AVG") & df$sequencing == "Metagenomic",])
four_mT = generate_cosm_area_plot(df[df$category %in% c("m78", "T0_AVG") & df$sequencing == "Metatranscriptomic",])
four_s = generate_cosm_area_plot(df[df$category %in% c("m78", "T0_AVG") & df$sequencing == "16S",])
  
fish_mG = generate_fish_barplot(df[df$category %in% c("Fish") & df$sequencing == "Metagenomic",])
fish_mT = generate_fish_barplot(df[df$category %in% c("Fish") & df$sequencing == "Metatranscriptomic",])
fish_s = generate_fish_barplot(df[df$category %in% c("Fish") & df$sample %in% c("F1", "F2", "F3", "F4", "F5", "F6") & df$sequencing == "16S",])

gridded = plot_grid(
  plot_grid(
    fish_s + theme(legend.position="none"),
          one_s + theme(legend.position="none"), two_s + theme(legend.position="none"), three_s + theme(legend.position="none"), four_s + theme(legend.position="none"),
    fish_mG + theme(legend.position="none"),
          one_mG + theme(legend.position="none"), two_mG + theme(legend.position="none"), three_mG + theme(legend.position="none"), four_mG + theme(legend.position="none"),
    fish_mT + theme(legend.position="none"),
          one_mT+ theme(legend.position="none"), two_mT + theme(legend.position="none"), three_mT+ theme(legend.position="none"), four_mT+ theme(legend.position="none"),
          nrow = 3, labels = c("A", "", "", "", "", "B", "", "", "", "", "C")), 
  
  get_legend(three_mT + theme(legend.box.margin = margin(0, 0, 0, 12))), rel_widths = c(3, .6))


svg("combined_read_taxonomy.svg", width = 10, height = 6)
gridded
dev.off()
