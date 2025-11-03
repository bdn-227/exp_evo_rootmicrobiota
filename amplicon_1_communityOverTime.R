
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
otu = read.delim("raw_data/amplicon_evo_otu.tsv") %>% column_to_rownames("strain")
design = read.delim("raw_data/amplicon_evo_design.tsv")



# subset the evolution experiment
idx = (design$compartment == "root" & design$sc != "mixed") | (design$compartment == "input")
designSubset = design[idx, ]
otu_subset = otu[, idx]



# calculate relative abundance
norm = apply(otu_subset, 2, function(x) x/sum(x))
merged_otu = norm %>%
    data.frame() %>%
    rownames_to_column("strain") %>%
            pivot_longer(
                cols = -strain,
                names_to = "sampleID",
                values_to = "RA"
            ) %>%
        left_join(design, by = "sampleID") %>%
        left_join(taxo_table, by = "strain") %>%
        data.frame()



# sort the strains by abundance
meanRA = merged_otu %>% group_by(plant, sc, strain, gen, family) %>% summarize(RA = mean(RA))
inputDf = merged_otu[merged_otu$compartment == "input", colnames(meanRA)]
meanRAdf = bind_rows(meanRA, inputDf)


orderDf = meanRAdf %>% group_by(family) %>% summarize(sum = sum(RA)) %>% arrange(sum)
orderDf2nd = orderDf
orderDf2nd$family = glue("{orderDf2nd$family}a")
orderDfOG = rbind(orderDf, orderDf2nd) %>% arrange(sum)
orderVec = orderDfOG$family


# make facette plots with all from one condition
scaling_factor=0.5
for (plant in c("col0", "gifu"))
{
    for (sc in c("at", "lj")) 
    {
        lines = unique(merged_otu[merged_otu$sc == sc & merged_otu$plant == plant, "line"])
        lines = as.character(1:20)[as.character(1:20) %in% lines]

        muller_list = list()
        for (line in lines) {
            
            # process the muller df
            muller_list[[line]] = convert_to_muller(line, merged_otu)
        }

        # concat together and plot
        muller_df = bind_rows(muller_list) %>% 
                        mutate(syncom = paste0("*", str_to_title(sc), "*-SC")) %>% 
                        mutate(line = paste0("population ", line)) %>% 
                        mutate(line = factor(line, levels=paste0("population ", lines)))

        # coloring of the strains
        color_family = color_family[color_family$Family != "unkown",]
        colors = color_family$Color
        names(colors) = color_family$Family

        if (!is.null(orderVec)) 
        {
            muller_df$Group_id = factor(muller_df$Group_id, levels = orderVec)
            idx = match(orderVec, names(colors))
            idx = idx[!is.na(idx)]
            colors = colors[idx]
        }

        # fix axis
        breaks = seq(1, max(muller_df$Generation), 2)
        limits = c(1, max(muller_df$Generation))
        xlab = breaks -1
        xlab[1] = "0"

        # text formatting
        levels = unique(str_to_title(muller_df$line))
        muller_df$line = factor(str_to_title(muller_df$line), levels = levels)

        p = ggplot(muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) + 
            geom_area() +
            guides(linetype = FALSE, color = FALSE, fill=FALSE) + 
            scale_y_continuous(expand = c(0, 0), labels = as.character(25 * (0:4)), name = "Relative abundance (%)") +
            scale_fill_manual(name = "Identity", values = colors) +
            scale_color_manual(values = colors) + 
            scale_x_continuous(expand = c(0, 0), labels = xlab, breaks = breaks, limits = limits) + 
            facet_grid(cols = vars(line)) +
            main_theme + 
            theme(text = element_text(family = "Helvetica"),
                  legend.position = "right",
                  panel.spacing = unit(1.5, "lines"),
                  plot.title = element_text(hjust = 0.5), 
                  panel.border = element_rect(color = "black", fill = NA, size = 1),
                  strip.background = element_rect(color = "white"),
                  strip.text = element_text(colour = 'black'),
                  strip.text.y = element_markdown(size = 14, margin = unit(c(0,0,0,20), "pt")),
                  plot.margin = margin(0,0.3,0,0.1, "cm")
                  ) + 
            xlab(expression("Plant cycle, "~italic(t)))
        
        save_fig(p, glue("figures/amplicon_mullerPlot_no_legend_{plant}_{sc}"), width=length(lines)*6*scaling_factor, height=6*scaling_factor)

    }
}

