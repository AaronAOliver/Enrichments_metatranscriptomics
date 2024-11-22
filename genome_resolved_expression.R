library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
options(scipen = 999)
setwd("Enrichments_metatranscriptomics")

# Import rsemession data, protein and MAG descriptions
rsem = read.csv("binned_expression_fixed.tab", sep = "\t")
cazyme_table = read.csv("cazyme.tab", sep = "\t")
sulfatase_table = read.csv("sulfatase.tab", sep = "\t")
camper_table = read.csv("camper.tab", sep = "\t")

bin_table = read.csv("protein_bins.tab", sep = "\t")
tax_table = read.csv("bin_taxonomy.tab", sep = "\t")
parent_bin = unique.data.frame(tax_table[c("taxonomy", "parent_bin")])

# Combine these into a single large table
rsem = merge(x = rsem, y = cazyme_table, by = "protein", all.x = TRUE)
rsem = merge(x = rsem, y = sulfatase_table, by = "protein", all.x = TRUE)
rsem = merge(x = rsem, y = camper_table, by = "protein", all.x = TRUE)

rsem = merge(x = rsem, y = bin_table, by = "protein", all.x = TRUE)
rsem = merge(x = rsem, y = tax_table, by = "bin", all.x = TRUE)

# Order the x-axis samples
mapping_rate = read.csv("mapping_rate.tab", sep = "\t")
sample_list = c("FISH_1", "FISH_2", "FISH_3", "FISH_4", "FISH_5", "FISH_6", "T0_A", "T0_B", "T0_C", "T0_D", "M30_T03", "M30_T08", "M30_T16", "M36_T03", "M36_T08", "M36_T16", "M75_T03", "M75_T08", "M75_T16", "M78_T03", "M78_T08", "M78_T16")
alginate_list = c("*PL6", "PL6_*", "*PL7", "*PL7_*", "*PL17*", "*PL32*", "*PL34*", "*PL36*", "*PL38*")
carrageenan_list = c("*S1_7*", "*S1_19*", "*GH16*", "*GH127*", "*GH2", "*GH2_*", "*GH42*")
fucoidan_list = c("*S1_15*",	"*S1_16*",	"*S1_17*",	"*GH29*",	"*GH95*"	,"*GH168*",	"*GH141*",	"*S1_25*",	"*GH107*")

# Normalize TPMs for mapping rates (only 20-60% of reads are MAG-resolved)
rsem_norm = rsem
for (samp in sample_list) {
  rsem_norm[[samp]] = rsem_norm[[samp]] * (filter(mapping_rate, sample==samp)$map_rate * 0.01)
}
rsem_norm = rsem_norm[rsem_norm$taxonomy != "Unclassified",]
rsem_norm$annotation = ifelse(!is.na(rsem_norm$sulfatase), rsem_norm$sulfatase, ifelse(!is.na(rsem_norm$camper), rsem_norm$camper, rsem_norm$CAZyme))

#### Set graphical parameters for different ggplots
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
                       xlab("Metatranscriptomic Sample"),
                       scale_y_continuous(labels = scales::percent))
  
# Pivot dataframe based on taxonomy
pivot_taxomomy = function(df, cutoff=0) {
  df = df %>% 
    pivot_longer(
      cols = sample_list, 
      names_to = "sample",
      values_to = "tpm"
    ) %>% aggregate(tpm ~ sample + taxonomy, FUN = sum)
  df$sample = factor(df$sample, levels = sample_list)
  df_specific = df_specific[df_specific$tpm_total > cutoff,]
  
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
  df_specific = df_specific[df_specific$tpm_total > cutoff,]
  
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
  

  

######## LOOK AT TOTAL EXPRESSION
total_expression_plot = ggplot(pivot_taxomomy(rsem_norm), 
                               aes(fill=taxonomy, y=tpm, x=sample)) +
  bar_graph_params + ylab("Genome-resolved TPM")

svg("figures/metaT_resolved_total.svg", width=10, height=6)
total_expression_plot
dev.off()



###### LOOK AT PROT-CLASS SPECIFIC EXPRESSION
total_cazy_expression_plot = ggplot(pivot_taxomomy(rsem_norm[is.na(rsem_norm$sulfatase) & !is.na(rsem_norm$CAZyme),]),
                                    aes(fill=taxonomy, y=tpm, x=sample)) + 
  bar_graph_params + ylab("Genome-resolved CAZyme TPM")

svg("figures/metaT_resolved_total_CAZyme.svg", width=10, height=6)
total_cazy_expression_plot 
dev.off()

total_sulf_expression_plot = ggplot(pivot_taxomomy(rsem_norm[!is.na(rsem_norm$sulfatase),]),
       aes(fill=taxonomy, y=tpm, x=sample)) + 
  bar_graph_params + ylab("Genome-resolved Sulfatase TPM")

svg("figures/metaT_resolved_total_sulf.svg", width=10, height=6)
total_sulf_expression_plot
dev.off()

total_camper_expression_plot = ggplot(pivot_taxomomy(rsem_norm[!is.na(rsem_norm$sulfatase),]),
                                    aes(fill=taxonomy, y=tpm, x=sample)) + 
  bar_graph_params + ylab("Genome-resolved Polyphenol-degrading TPM")

svg("figures/metaT_resolved_total_camper.svg", width=10, height=6)
total_camper_expression_plot
dev.off()

####### DEREPLICATED EXPRESSION AVERAGES - SULFATASE

rsem_total_derep = rsem_norm %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_total"
  ) %>% aggregate(tpm_total ~ sample + parent_bin, FUN = sum)

rsem_sulf_derep = rsem_norm[!is.na(rsem_norm$sulfatase), ] %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_relevant"
  ) %>% aggregate(tpm_relevant ~ sample + parent_bin, FUN = sum)

rsem_total_derep = merge(x = rsem_total_derep, y = parent_bin, by = "parent_bin", all.x = TRUE)
rsem_sulf_derep = merge(x = rsem_sulf_derep, y = rsem_total_derep, by = c("sample","parent_bin"), all.x = TRUE)

rsem_sulf_derep = rsem_sulf_derep[rsem_sulf_derep$tpm_total > 5000,]
rsem_sulf_derep$pct_relevant = (rsem_sulf_derep$tpm_relevant + 1) / (rsem_sulf_derep$tpm_total + 1)

svg("figures/metaT_resolved_sulf_pct.svg", width=10, height=6)
ggplot(rsem_sulf_derep, aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
  geom_beeswarm(cex = 0.5, size = 2) + theme_bw() + 
  scale_color_manual(name = "MAG Taxonomy", values = taxonomy_coloring) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")) + 
  ylab("Sulfatases as % of MAG-specific Expression") + xlab("Metatranscriptomic Sample") +
  scale_y_continuous(labels = scales::percent) + scale_x_discrete(drop = FALSE)
dev.off()

######## DEREPLICATED EXPRESSION AVERAGES - CAZYME

rsem_total_derep = rsem_norm %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_total"
  ) %>% aggregate(tpm_total ~ sample + parent_bin, FUN = sum)

rsem_cazy_derep = rsem_norm[is.na(rsem_norm$sulfatase) & !is.na(rsem_norm$CAZyme), ] %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_relevant"
  ) %>% aggregate(tpm_relevant ~ sample + parent_bin, FUN = sum)

rsem_total_derep = merge(x = rsem_total_derep, y = parent_bin, by = "parent_bin", all.x = TRUE)
rsem_cazy_derep = merge(x = rsem_cazy_derep, y = rsem_total_derep, by = c("sample","parent_bin"), all.x = TRUE)

rsem_cazy_derep = rsem_cazy_derep[rsem_cazy_derep$tpm_total > 5000,]
rsem_cazy_derep$pct_relevant = (rsem_cazy_derep$tpm_relevant + 1) / (rsem_cazy_derep$tpm_total + 1)

svg("figures/metaT_resolved_CAZyme_pct.svg", width=10, height=6)
ggplot(rsem_cazy_derep, aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
  geom_beeswarm(cex = 0.5, size = 2) + theme_bw() + 
  scale_color_manual(name = "MAG Taxonomy", values = taxonomy_coloring) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")) + 
  ylab("CAZymes as % of MAG-specific Expression") + xlab("Metatranscriptomic Sample") +
  scale_y_continuous(labels = scales::percent)
dev.off()


###### HIGHEST INDIVIDUAL GENE EXPRESSION - ALGINATE

rsem_alginate_derep = rsem_norm[Reduce("|", lapply(alginate_list, grepl, rsem_norm$CAZyme)), ] %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_total"
  ) %>% aggregate(tpm_total ~ sample + parent_bin + CAZyme, FUN = sum)
rsem_alginate_derep = merge(x = rsem_alginate_derep, y = parent_bin, by = "parent_bin", all.x = TRUE)
rsem_alginate_derep = rsem_alginate_derep[rsem_alginate_derep$tpm_total > 5,]

svg("figures/metaT_resolved_total_alginate.svg", width=10, height=6)
ggplot(rsem_alginate_derep, aes(x = sample, y = tpm_total, color=taxonomy, group=parent_bin)) +
  geom_beeswarm(cex = 0.5, size = 2) + theme_bw() + 
  scale_color_manual(name = "MAG Taxonomy", values = taxonomy_coloring) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")) + 
  ylab("Alginate Enzyme TPM") + xlab("Metatranscriptomic Sample")
dev.off()

rsem_total_derep = rsem_norm %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_total"
  ) %>% aggregate(tpm_total ~ sample + parent_bin, FUN = sum)
rsem_total_derep = merge(x = rsem_total_derep, y = parent_bin, by = "parent_bin", all.x = TRUE)

rsem_alginate_pct_derep = rsem_norm[Reduce("|", lapply(alginate_list, grepl, rsem_norm$annotation)), ]  %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_relevant"
  ) %>% aggregate(tpm_relevant ~ sample + parent_bin, FUN = sum)
rsem_alginate_pct_derep$sample = factor(rsem_alginate_pct_derep$sample, levels = sample_list)

rsem_alginate_pct_derep = merge(x = rsem_alginate_pct_derep, y = rsem_total_derep, by = c("sample","parent_bin"), all.x = TRUE)

rsem_alginate_pct_derep = rsem_alginate_pct_derep[rsem_alginate_pct_derep$tpm_total > 5000,]
rsem_alginate_pct_derep$pct_relevant = (rsem_alginate_pct_derep$tpm_relevant + 1) / (rsem_alginate_pct_derep$tpm_total + 1)

svg("figures/metaT_resolved_alginate_pct.svg", width=10, height=6)
ggplot(rsem_alginate_pct_derep, aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
  geom_beeswarm(cex = 0.5, size = 2) + theme_bw() + scale_x_discrete(drop = FALSE) +
  scale_color_manual(name = "MAG Taxonomy", values = taxonomy_coloring) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")) + 
  ylab("Alginate-degradation as % of MAG-specific Expression") + xlab("Metatranscriptomic Sample") +
  scale_y_continuous(labels = scales::percent)
dev.off()

######

rsem_carra_derep = rsem_norm[Reduce("|", lapply(carrageenan_list, grepl, rsem_norm$annotation)), ] %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_total"
  ) %>% aggregate(tpm_total ~ sample + parent_bin + annotation, FUN = sum)
rsem_carra_derep = merge(x = rsem_carra_derep, y = parent_bin, by = "parent_bin", all.x = TRUE)
rsem_carra_derep = rsem_carra_derep[rsem_carra_derep$tpm_total > 10,]

svg("figures/metaT_resolved_total_carra.svg", width=10, height=6)
ggplot(rsem_carra_derep, aes(x = sample, y = tpm_total, color=taxonomy)) +
  geom_beeswarm(cex = 0.5, size = 2) + theme_bw() + 
  scale_color_manual(name = "MAG Taxonomy", values = taxonomy_coloring) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")) + 
  ylab("Carrageenan Enzyme TPM") + xlab("Metatranscriptomic Sample")
dev.off()

pct_carrageenan_expression_plot = ggplot(pivot_parent_bin(rsem_norm, rsem_norm[Reduce("|", lapply(carrageenan_list, grepl, rsem_norm$annotation)), ], 5000),
                   aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
  beeswarm_params + ylab("Genome-resolved Sulfatase TPM")

svg("figures/metaT_resolved_carra_pct.svg", width=10, height=6)
pct_carrageenan_expression_plot
dev.off()



#####


rsem_fucoidan_derep = rsem_norm[Reduce("|", lapply(fucoidan_list, grepl, rsem_norm$annotation)), ] %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_total"
  ) %>% aggregate(tpm_total ~ sample + parent_bin + annotation, FUN = sum)
rsem_fucoidan_derep = merge(x = rsem_fucoidan_derep, y = parent_bin, by = "parent_bin", all.x = TRUE)
rsem_fucoidan_derep = rsem_fucoidan_derep[rsem_fucoidan_derep$tpm_total > 10,]

svg("figures/metaT_resolved_total_fcsp.svg", width=10, height=6)
ggplot(rsem_fucoidan_derep, aes(x = sample, y = tpm_total, color=taxonomy)) +
  geom_beeswarm(cex = 0.5, size = 2) + theme_bw() + 
  scale_color_manual(name = "MAG Taxonomy", values = taxonomy_coloring) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")) + 
  ylab("Fucoidan Enzyme TPM") + xlab("Metatranscriptomic Sample")
dev.off()

rsem_total_derep = rsem_norm %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_total"
  ) %>% aggregate(tpm_total ~ sample + parent_bin, FUN = sum)

rsem_fucoidan_pct_derep = rsem_norm[Reduce("|", lapply(fucoidan_list, grepl, rsem_norm$annotation)), ]  %>% 
  pivot_longer(
    cols = sample_list, 
    names_to = "sample",
    values_to = "tpm_relevant"
  ) %>% aggregate(tpm_relevant ~ sample + parent_bin, FUN = sum)

rsem_total_derep = merge(x = rsem_total_derep, y = parent_bin, by = "parent_bin", all.x = TRUE)
rsem_fucoidan_pct_derep = merge(x = rsem_fucoidan_pct_derep, y = rsem_total_derep, by = c("sample","parent_bin"), all.x = TRUE)

rsem_fucoidan_pct_derep = rsem_fucoidan_pct_derep[rsem_fucoidan_pct_derep$tpm_total > 5000,]
rsem_fucoidan_pct_derep$pct_relevant = (rsem_fucoidan_pct_derep$tpm_relevant + 1) / (rsem_fucoidan_pct_derep$tpm_total + 1)

svg("figures/metaT_resolved_fcsp_pct.svg", width=10, height=6)
ggplot(rsem_fucoidan_pct_derep, aes(x = sample, y = pct_relevant, color=taxonomy, group=parent_bin)) +
  geom_beeswarm(cex = 0.5, size = 2) + theme_bw() +
  scale_color_manual(name = "MAG Taxonomy", values = taxonomy_coloring) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(color="black")) + 
  ylab("Fucoidan-degradation as % of MAG-specific Expression") + xlab("Metatranscriptomic Sample") +
  scale_y_continuous(labels = scales::percent)
dev.off()

############

assembled_rsem = read.csv("assembled_expression.tab", sep = "\t")  %>% 
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
