###Summarised workflow of R analysis
#Paper: "From Forest to Farm: The Impact of a Broad Spectrum of Lifestyles on the Porcine Gut Microbiota", Comer et al. 2025

##Inputs (also available via Zeonodo)
#ps - phyloseq object of all data 
#ps.cor - phyloseq object after batch-effect correction (sPLSDA-batch method: )
#ps.dom - phyloseq object of all non-wild animals (wild boar removed)
#ps.cor.dom - phyloseq object of all non-wild animals (after batch-effect correction)
#genus.tab - genus table (relative abundances)
#maastab - genus table (raw reads, used for MaAsLin2)
#metadata - the metadata input used for MaAsLin2


library(phyloseq)          
library(vegan)             
library(Maaslin2)          
library(microbiomeMarker)  
library(dplyr)             
library(tidyr)             
library(knitr)           
library(ggplot2)           
library(ggpubr)            
library(ggExtra)          
library(dominanceanalysis) 
library(randomForest)    
library(scales)           



###########################
##Diversity area modelling 

sample_sizes <- seq(5, ncol(genus.tab), by = 10)
diversity_results <- data.frame(Samples = integer(), Richness = integer())
set.seed(123)
for (n in sample_sizes) {
  richness_vals <- replicate(50, {
    selected_samples <- sample(colnames(genus.tab), n)
    subsample_data <- genus.tab[, selected_samples]
    sum(rowSums(subsample_data) > 0)  # count genera present
  })
  diversity_results <- rbind(diversity_results, data.frame(
    Samples = n,
    Richness = mean(richness_vals)
  ))
}
#Fit a power-law model: S = c*A^z
sar_model <- nls(Richness ~ c * Samples^z, data = diversity_results, start = list(c = 100, z = 0.3))
summary(sar_model)
plot(Richness ~ Samples, data = diversity_results, pch = 16, main = "Diversity–Area Curve")+
  lines(diversity_results$Samples, predict(sar_model), col = "blue", lwd = 2)








###########################
##Calculating the core microbiota
#Set threshold (in this case, 95%)
threshold <- 0.95
sample_count <- ncol(genus.tab)
cutoff <- threshold * sample_count

#Count non-zero abundances per genus
present_in_samples <- rowSums(genus.tab > 0)

#Filter genera present in more than 95% of samples
genera_95 <- rownames(genus.tab)[present_in_samples > cutoff]
genera_95

core_abundance <- genus.tab[genera_95, ]
core_sum_per_sample <- colSums(core_abundance)

#Get the median core abundance across all samples
median_core_abundance <- median(core_sum_per_sample)
median_core_abundance




###########################
##Alpha-diversity 
#Extract the OTU table
otu_data <- otu_table(ps.dom)

#Calculate Pielou's Evenness for each sample
Pielou_evenness <- apply(otu_data, 1, function(x) {
  # Shannon diversity (H')
  H_prime <- diversity(x, index = "shannon")
  # Richness (S) (number of non-zero species)
  S <- length(x[x > 0])
  # Pielou's Evenness (J')
  J_prime <- H_prime / log(S)
  return(J_prime)
})

#Add the Pielou's Evenness values as a new column to the sample data
sample_data(ps.dom)$Pielou_evenness <- Pielou_evenness

#Extract the alpha-diversity data (Simpson, Shannon, Chao1) and Pielou's Evenness
#Create a custom data frame that includes the diversity measures
alpha_diversity_data <- data.frame(
  Group = rep(sample_data(ps.dom)$Age_category, times = 4),  # Repeat for 4 measures
  Measure = rep(c("Simpson", "Shannon", "Chao1", "Pielou Evenness"), each = length(sample_data(ps.dom)$Age_category)),
  Value = c(
    estimate_richness(ps.dom, measures = "Simpson")[,1],  # Simpson
    estimate_richness(ps.dom, measures = "Shannon")[,1],  # Shannon
    estimate_richness(ps.dom, measures = "Chao1")[,1],    # Chao1
    sample_data(ps.dom)$Pielou_evenness                  # Pielou's Evenness
  )
)

group_colours <- c(
  "Lactation" = "#e66500",   
  "Nursery"   = "#56B4E9",   
  "Growing"   = "#0bb071",   
  "Finishing" = "#F0E442",   
  "Mature"    = "#ca3878"    
)
# Plot
ggplot(alpha_diversity_data, aes(x = Group, y = Value, fill = Group)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.2) +
  geom_boxplot(size = 0.8, staplewidth = 0.2) +
  scale_fill_manual(values = group_colours) +
  stat_compare_means(method = "kruskal.test", size = 5) +
  facet_wrap(~Measure, scales = "free", ncol = 4) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "#f8faf2"),
    #axis.text.x = element_blank(),
    strip.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 20.2),
    legend.text = element_text(size = 20),
    axis.title.y = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.key.size = unit(2, "line")
  ) +
  labs(y = "Alpha Diversity Measure")








###########################
##Beta-diversity 

#Transform data and choose metric (in this case Bray-Curtis)
ps_log <- transform_sample_counts(ps.cor, function(x) log(1+x))
pcoa<-ordinate(ps_log, "MDS", "bray" )

# Define color and shape palettes
group_colours <- c(
  "Lactation" = "#e66500",   
  "Nursery"   = "#56B4E9",   
  "Growing"   = "#0bb071",   
  "Finishing" = "#F0E442",   
  "Mature"    = "#ca3878"    
)

shapes <- c(16, 17, 15, 18, 16)  
fill_palette <- group_colours  

#Base ordination plot
pcoa_plot <- plot_ordination(ps_log, pcoa, type = "sample", color = "Age_category", shape = "Age_category") +
  geom_point(size = 5) +  
  stat_ellipse(aes(group = Age_category, fill = factor(Age_category)), 
               type = "norm", level = 0.9, lty = 0,
               geom = "polygon", alpha = 0.05, show.legend = FALSE) + 
  scale_color_manual(values = group_colours) +
  scale_shape_manual(values = shapes) +
  scale_fill_manual(values = fill_palette) +
  ggtitle("Age category") +
  labs(x = "PCo1 (8.6%)", y = "PCo2 (6.7%)", 
       color = "Age category", shape = "Age category") +
  theme_light() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    panel.border = element_rect(color = "black", size = 1.4, fill = NA),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  )

#Add density plots
pcoa_plot_bray <- ggMarginal(
  pcoa_plot,
  type = "density",
  margins = "both",
  groupColour = TRUE,
  groupFill = TRUE
)

#Final plot
pcoa_plot_bray





###########################
##LEfSe analysis 
table(sample_data(ps.dom)$Age_category)
lef_out<-run_lefse(ps.dom, group = "Age_category", norm = "CPM",
                   taxa_rank = "Genus", multigrp_strat = TRUE,
                   
                   kw_cutoff = 0.05, lda_cutoff = 2.5)

lef_out@marker_table
dat<-marker_table(lef_out) %>% data.frame() %>% select(1:5)
print(dat)
dat = dat %>% mutate(p.experiment =p.adjust(pvalue, method="BH"))
dat %>% kable(align = "c")

#write.csv(dat, "c:/Users/Location/Results_LEfSe.csv")
plot_ef_bar(lef_out)









###########################
##MaAsLin2

fit_data <- Maaslin2(
  input_data = maastab,
  input_metadata = metadata,
  output = "maaslin2_output",
  fixed_effects = c("Crude.fibre", "Slatted.floor", "Bedding_present", "Breed_type", 
                    "Soil.present", "Outdoor_access"),
  normalization = "TSS",      #Or CLR, 
  transform = "LOG",          
  analysis_method = "LM",     #Linear model (for continuous + categorical)
  standardize = TRUE
)

results <- read.table("maaslin2_output/significant_results.tsv", sep = "\t", header = TRUE)
head(results)





###########################
##Evaluating the effects of different factors 

###PERMANOVA and omega-squared effect size
#Extract from phyloseq the tables
seqtab_log<-as(otu_table(ps_log), "matrix")
taxonomy_log<-as(tax_table(ps_log), "matrix")
design_log<-as(sample_data(ps_log), "data.frame")
#dissimilarity matrix
dist.bray<-vegdist(seqtab_log, method = "bray")
#define variable
groups <- as.factor(design_log$Age_category)
levels(groups)
#Basic PERMANOVA test
#adonis2(dist.bray ~ groups, design_log, permutations = 9999)

#Omega-squared
per<-adonis2(dist.bray ~ groups, design_log, permutations = 9999)
adonis_OmegaSq(per, partial = TRUE)









###
#Dominance analysis 
Bray-Curtis distance and PCoA
pcoa <- ordinate(ps_log, method = "PCoA", distance = "bray")
pc1 <- pcoa$vectors[,1]  # Extract PC1 values
#Combine metadata with PC1
design_log$PC1 <- pc1
#Run linear model
lm_model <- lm(PC1 ~ Age_category + Country + Species + 
                 Breed_type + Husbandry.intensity + Outdoor_access + 
                 Soil.present + Bedding.material + Bedding_present + 
                 Straw_present + Slatted.floor + Purely_commercial_diet +
                 Seq_run + CH_Enterotype + DMM_cluster + HDBSCAN, data = design_log)
#Run dominance analysis
da <- dominanceAnalysis(lm_model)
summary(da)
plot(da)




###
#Random forest 
#Collapse to genus level
ps_genus <- tax_glom(ps.dom, taxrank = "Genus")
#Filter rare genera (present in >5% samples)
ps_genus <- filter_taxa(ps_genus, function(x) sum(x > 0) > 0.05 * length(x), TRUE)
#Convert counts to relative abundance
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
#Replace taxa names with genus names
genus_names <- as.character(tax_table(ps_genus_rel)[, "Genus"])
genus_names[is.na(genus_names)] <- "Unclassified_Genus"
genus_names <- make.unique(genus_names)
taxa_names(ps_genus_rel) <- genus_names
#Extract OTU table as a data frame with samples as rows, genera as columns
otu_mat <- otu_table(ps_genus_rel)
if (taxa_are_rows(otu_mat)) {
  otu_mat <- t(otu_mat)
}
otu_data <- as.data.frame(otu_mat)
otu_data$sample <- rownames(otu_data)
#Extract metadata
meta_data <- as(sample_data(ps_genus_rel), "data.frame")
meta_data$sample <- rownames(meta_data)
#Merge on sample (now samples are rows in both data frames)
merged <- merge(meta_data, otu_data, by = "sample")
#Confirm Age_category is factor
merged$Age_category <- as.factor(merged$Age_category)
# Prepare data for random forest:
# - response variable: Age_category
# - predictor variables: genera abundance columns (all columns after metadata columns)
# Note: adjust 'metadata_end' if metadata columns count is different
metadata_cols <- colnames(meta_data)
metadata_end <- max(match(metadata_cols, colnames(merged)))
genus_start <- metadata_end + 1
rf_input <- merged[, c("Age_category", colnames(merged)[genus_start:ncol(merged)])]
# Make column names syntactically valid
colnames(rf_input) <- make.names(colnames(rf_input))
# Remove rows with missing values
rf_input <- na.omit(rf_input)
# Train random forest model
set.seed(42)
rf_model <- randomForest(Age_category ~ ., data = rf_input, importance = TRUE, ntree = 1000)
# Check results
print(rf_model)
varImpPlot(rf_model, n.var = 20, main = "Top Predictive Genera")



###########################


