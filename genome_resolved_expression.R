## Set up environment
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
options(scipen = 999)
setwd("Enrichments_metatranscriptomics")

# Import Expression data, protein and MAG descriptions
mapping_rate = read.csv("mapping_rate.tab", sep = "\t")
tax_table = read.csv("bin_taxonomy.tab", sep = "\t")
parent_bin = unique.data.frame(tax_table[c("taxonomy", "parent_bin")])

rsem = left_join(read.csv("binned_expression_fixed.tab", sep = "\t"),
                 read.csv("cazyme.tab", sep = "\t"), by= "protein") %>%
       left_join(., read.csv("sulfatase.tab", sep = "\t"), by = "protein") %>%
       left_join(., read.csv("camper.tab", sep = "\t"), by = "protein") %>%
       left_join(., read.csv("protein_bins.tab", sep = "\t"), by = "protein") %>%
       left_join(., tax_table, by = "bin")
  
# Order the x-axis samples and categories of protein regex
sample_list = c("FISH_1", "FISH_2", "FISH_3", "FISH_4", "FISH_5", "FISH_6", "T0_A", "T0_B", "T0_C", "T0_D", "M30_T03", "M30_T08", "M30_T16", "M36_T03", "M36_T08", "M36_T16", "M75_T03", "M75_T08", "M75_T16", "M78_T03", "M78_T08", "M78_T16")
alginate_list = c("*PL6", "PL6_*", "*PL7", "*PL7_*", "*PL17*", "*PL32*", "*PL34*", "*PL36*", "*PL38*")
carrageenan_list = c("*S1_7*", "*S1_19*", "*GH16*", "*GH127*", "*GH2", "*GH2_*", "*GH42*")
fucoidan_list = c("*S1_15*",	"*S1_16*",	"*S1_17*",	"*GH29*",	"*GH95*"	,"*GH168*",	"*GH141*",	"*S1_25*",	"*GH107*")

# Normalize TPMs for mapping rates (only 20-60% of reads are MAG-resolved)
rsem_norm = rsem
for (samp in sample_list) {
  rsem_norm[[samp]] = rsem_norm[[samp]] * (filter(mapping_rate, sample==samp)$map_rate * 0.01)
}

# Drop poorly classified MAGs and standardize annotation column
rsem_norm = rsem_norm[rsem_norm$taxonomy != "Unclassified",]
rsem_norm$annotation = ifelse(!is.na(rsem_norm$sulfatase), rsem_norm$sulfatase, ifelse(!is.na(rsem_norm$camper), rsem_norm$camper, rsem_norm$CAZyme))

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
  
# Functions to pivot dataframe or dataframe pair
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

## Total Expression Figure
total_expression_plot = ggplot(pivot_taxomomy(rsem_norm), 
                               aes(fill=taxonomy, y=tpm, x=sample)) +
  bar_graph_params + ylab("Genome-resolved TPM")

svg("figures/metaT_resolved_total.svg", width=10, height=6)
total_expression_plot
dev.off()

## Bulk CAZyme Expression
{
total_cazy_expression_plot = ggplot(pivot_taxomomy(rsem_norm[is.na(rsem_norm$sulfatase) & !is.na(rsem_norm$CAZyme),]),
                                    aes(fill=taxonomy, y=tpm, x=sample)) + 
  bar_graph_params + ylab("Genome-resolved CAZyme TPM")

svg("figures/metaT_resolved_total_CAZyme.svg", width=10, height=6)
total_cazy_expression_plot 
dev.off()

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
total_sulf_expression_plot = ggplot(pivot_taxomomy(rsem_norm[!is.na(rsem_norm$sulfatase),]),
       aes(fill=taxonomy, y=tpm, x=sample)) + 
  bar_graph_params + ylab("Genome-resolved Sulfatase TPM")

svg("figures/metaT_resolved_total_sulf.svg", width=10, height=6)
total_sulf_expression_plot
dev.off()

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
total_camper_expression_plot = ggplot(pivot_taxomomy(rsem_norm[!is.na(rsem_norm$sulfatase),]),
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

##### Expression based on metaSPAdes assembled metatranscripts
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
