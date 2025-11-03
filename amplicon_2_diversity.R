# script by Niklas Kiel
# community profiles over time
options(warn=-1)

# load packages
library(ggplot2)
library(tidyverse)
library(reshape2)
library(psych)
library(tidyr)
library(glue)
library(ggtext)


# load functions
source("scripts/amplicon_analysis_functions.R")
source("scripts/plotting_functions.R")


# read the otu table and design
otu = read.delim("tables/amplicon_evo_otu.tsv")
design = read.delim("tables/amplicon_evo_design.tsv")


# determine variance explained
# get subsets
idx = (design$compartment == "root" & design$sc != "mixed")
designSubset = design[idx, ]
otuSubset = otu[, idx]

# get otu table and design
normalized_otu = t(otuSubset) / colSums(otuSubset)
dist = vegdist(normalized_otu)
rowidx = match(rownames(normalized_otu), design$SampleID)
d = design[rowidx, c("SampleID", "gen", "plant", "sc")]

# use cap-scale to identify the total variance
cap_plant_partial = capscale(dist ~ plant + Condition(sc) + Condition(gen), data = d)
pval = anova.cca(cap_plant_partial, permutations=9999)["Pr(>F)"]
cons_eig = cap_plant_partial$CCA$eig
uncons_eig = cap_plant_partial$CA$eig
sum_cons = sum(cons_eig, na.rm = TRUE)
sum_uncons = sum(uncons_eig, na.rm = TRUE)
total = sum_cons + sum_uncons
pct_cons = sum_cons / total * 100
cat("Plant (partial) explains", round(pct_cons, 2), "% of total inertia\n")




# within each host species
taxo_table = taxo_table[taxo_table$strain %in% rownames(otuSubset), ]
df_ls = list()
for (family in unique(taxo_table$family))
{
    strains = taxo_table[taxo_table$family == family, ]$strain
    df = data.frame(colSums(otuSubset[strains, ]))
    colnames(df) = family
    df_ls[[family]] = df
}

# perform the analysis
otu_ = t(bind_cols(df_ls))
normalized_otu = t(otu_) / colSums(otu_)
dist = vegdist(normalized_otu)
rowidx = match(rownames(normalized_otu), design$SampleID)
d = design[rowidx, c("SampleID", "gen", "plant", "sc", "line")]
adonis_res = adonis2(dist ~ gen, data=d, permutations=9999, by = "margin")
print(adonis_res)




# calculate percent incerase for flavobacteria
strains = taxo_table[taxo_table$strain == "Root935", ]$strain
idx = design$compartment == "root" & design$sc != "mixed"
designSubset = design[idx, ]
otuSubset = otu[, idx]
norm = apply(otu, 2, function(x) x/sum(x))
melted = melt(norm)
colnames(melted) = c("strain", "SampleID", "RA")
mergedOtuTab = melted %>% left_join(design) %>% filter(strain %in% strains)



percent_increase <- function(old, new) {
  if (old == 0) {
    stop("Old value cannot be zero when calculating percent increase.")
  }
  ((new - old) / old) * 100
}

# calculcate percent increase
round(percent_increase(mean(mergedOtuTab[mergedOtuTab["gen"]==1, "RA"]), mean(mergedOtuTab[mergedOtuTab["gen"]==16, "RA"])), 2)

