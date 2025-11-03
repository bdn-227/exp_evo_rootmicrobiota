
# packages
library(ggmuller)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(vegan)
library(ggplot2)





variability_table <- function(cca){
  
  chi <- c(cca$tot.chi,
           cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi/chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained")
  return(variability_table)
  
}

cap_var_props <- function(cca){
  
  eig_tot <- sum(cca$CCA$eig)
  var_propdf <- cca$CCA$eig/eig_tot
  return(var_propdf)
}

pca_var_props <- function(cca){
  
  eig_tot <- sum(cca$CA$eig)
  var_propdf <- cca$CA$eig/eig_tot
  return(var_propdf)
}

cca_ci <- function(cca, permutations=5000){
  
  var_tbl <- variability_table(cca)
  p <- permutest(cca, permutations=permutations)
  ci <- quantile(p$F.perm, c(.05,.95))*p$chi[1]/var_tbl["total", "inertia"]
  return(ci)
  
}




# convert abundance data to muller df
convert_to_muller = function(line, merged_otu)
{
    line = as.character(line)
    merged_muller_df = data.frame(merged_otu)
    sc = unique(merged_muller_df[merged_muller_df$line == line, "sc"])
    host_plant = unique(merged_muller_df[merged_muller_df$line == line, "plant"])
    idx = (merged_muller_df$line == line | (merged_muller_df$compartment == "input")) & merged_muller_df$sc == sc
    mullerDf = merged_muller_df[idx, ]
    mullerDf$gen = as.numeric(mullerDf$gen)

    # remove non-syncom strains
    if (sc == "at") keepStr = "Root"
    if (sc == "lj") keepStr = "Lj"
    mullerDf = mullerDf[str_starts(mullerDf$strain, keepStr) | str_starts(mullerDf$strain, "LjNodule"), ]

    # shift generations 
    mullerDf$gen = mullerDf$gen + 1


    edges = data.frame(Parent = "origin", Identity= unique(mullerDf$family))
    origin_df = data.frame(SampleID = NA, 
                            sample_name=NA,
                            strain = NA, 
                            RA = 1, 
                            sc = NA, 
                            plant = NA, 
                            compartment = NA, 
                            line = NA, 
                            gen = 0, 
                            family = "origin")

    # add missing columns back
    missingCols = colnames(origin_df)[!colnames(origin_df) %in% colnames(mullerDf)]
    mullerDf[missingCols] = NA
    mullerDf = mullerDf[, colnames(origin_df)]

    pop_df = rbind(mullerDf, origin_df)
    pop_df = pop_df[, c("gen", "family", "RA")]
    colnames(pop_df) = colnames(example_pop_df)
    muller_df = get_Muller_df(edges, pop_df %>% group_by(Generation, Identity) %>% summarize(Population = sum(Population)))

    # subset only strains
    muller_df_strains = muller_df %>% 
                            dplyr::filter(Identity != "origin") %>% 
                            mutate(line = !!line)
    return(muller_df_strains)
}



# function to plot profiles over time
plot_muller = function(line, merged_muller_df, orderVec = NULL) {

    muller_df_strains = convert_to_muller(line, merged_muller_df)

    # coloring of the strains
    color_family = color_family[color_family$Family != "unkown",]
    colors = color_family$Color
    names(colors) = color_family$Family

    # gather meta-information
    title_string = paste0("line ", line)

    if (!is.null(orderVec)) {
      muller_df_strains$Group_id = factor(muller_df_strains$Group_id, levels = orderVec)
      idx = match(orderVec, names(colors))
      idx = idx[!is.na(idx)]
      colors = colors[idx]
    }

    # fix axis
    breaks = seq(1, max(muller_df_strains$Generation), 2)
    limits = c(1, max(muller_df_strains$Generation))
    xlab = breaks -1
    xlab[1] = "In"

    p = ggplot(muller_df_strains, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) + 
          geom_area() +
          theme(legend.position = "right") +
          guides(linetype = FALSE, color = FALSE) + 
          scale_y_continuous(expand = c(0, 0), labels = 25 * (0:4), name = "Relative abundance (%)") +
          scale_fill_manual(name = "Identity", values = colors) +
          scale_color_manual(values = colors) + 
          facette_theme + 
          scale_x_continuous(expand = c(0, 0), labels = xlab, breaks = breaks, limits = limits)
    
    return(p)
}


