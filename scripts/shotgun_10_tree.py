

# load module
from Bio import Phylo
import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns
from matplotlib import ticker
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
import matplotlib.colors as mcolors
import numpy as np

# load custom stuff
from shotgun_analysis_functions import load_colors, color_dict_strain
plt.rcParams.update({
    "font.family": "Helvetica"
})

# strings
title_dict = { 
            'at(gifu)':  "$At$-SC$^{Lj}$", 
            'at(col0)':  "$At$-SC$^{At}$",
            'lj(col0)':  "$Lj$-SC$^{At}$",
            'lj(gifu)':  "$Lj$-SC$^{Lj}$", 
            }


# load data
evolved_table = pl.read_csv("mapping/evolved_strains.tsv", 
                            separator="\t")
color_dict = color_dict_strain()
colors_all = load_colors()
tree_colors = (pl.read_csv("results/tree_colors.tsv", separator="\t")
                .with_columns(pl.col("color").alias("edgecolor"))
                .with_columns(pl.col("color").str.replace("black", "lightgrey"))
                )

# LOAD THE GENOMIC SIZES
genome_sizes = (pl.read_csv('raw_data/coverage_strainwise.tsv', separator="\t")
                . with_columns( ("contig_" + pl.col("#rname").str.split_exact("_", n=2).struct[-1]).alias("contig") )
                .select("strain", "contig", "endpos")
                .unique()
                .sort("strain", "contig")
                .with_columns(pl.col("endpos")/1e6)
                )

colors = {"$A. thaliana$ (Col-0)": "red",
          "$L. japonicus$ (Gifu)": "blue",
          "Trajectory": "grey",
          "Unspecific": "black"
          }




# ~~~~~~~~~~~~~~~~ PLOTTING ALL HITS ~~~~~~~~~~~~~~~~ #
variants = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv", 
                            separator="\t",
                            schema_overrides={"line": str}
                            ) 
                            .filter(pl.col("line").is_in(["1","2","3","4","5",
                                                          "6.2","7.2","8.2","9","10",
                                                          "11","12","13","14","15",
                                                          "16","17","18","19","20"]))
                            .select("ID", "syncom", "plant", "line", "AF") 
                            .group_by("ID", "syncom", "plant").agg(pl.col("AF").max())
                            .rename({"plant": "enriched"})
                            .with_columns(pl.col("ID").str.split_exact("_", n=5).struct.rename_fields(["strain", "contig", "contig_n", "loc", "ref", "mut"]).alias("fields"))
                            .unnest("fields")
                            .with_columns((pl.col("contig") + "_" + pl.col("contig_n")).alias("contig"))
                            .with_columns(pl.col("loc").cast(float)/1e6)
                            .with_columns(pl.lit(True).alias("significant"))
                            )
# ~~~~~~~~~~~~~~~~ PLOTTING ALL HITS ~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~ PLOTTING ALL HITS + PACBIO ~~~~~~~~~~~~~~~~ #
variants_pacbio = (pl.read_csv('results/shotgun_pacbio_all.tsv', separator="\t", schema_overrides={"line": str})
                        .with_columns(pl.col("ID").str.split_exact("_", n=5).struct.rename_fields(["strain", "contig", "contig_n", "loc", "ref", "mut"]).alias("fields"))
                        .unnest("fields")
                        .with_columns((pl.col("contig") + "_" + pl.col("contig_n")).alias("contig"))
                        .with_columns(pl.col("loc").cast(float)/1e6)
                        .with_columns( pl.when(pl.col("strain").str.starts_with("LjRoot")).then(pl.lit("lj"))
                                        .when(pl.col("strain").str.starts_with("Root")).then(pl.lit("at")) 
                                        .otherwise("syncom").alias("syncom")
                                        )
                        .with_columns( (pl.col("strain") + "-" + pl.col("plant") + "-" + pl.col("syncom")).alias("strain-idx") )
                        .sort("plant", "syncom", "line", "strain")
                        .filter(pl.col("af_shotgun").is_not_null())
                        .filter(pl.col("af_shotgun") >= 0.95)
                        )


# define palette
palette = {"col0": "red", "gifu": "blue"}
markers = {"native": "s", "non-native": "D"}


for syncom in ["at", "lj"]: 

    # subset the evolved strains
    evolved_sc = evolved_table.filter(pl.col("sc_origin") == syncom).sort("condition", descending=True)
    n_syncoms = evolved_sc["syncom"].unique().len()

    # load tree file
    tree = Phylo.read(f"mapping/species_tree_{syncom}.txt", "newick")
    tree.rooted = True

    # create the figure
    fig, axes = plt.subplots(figsize=(12,5), 
                             ncols=3,)

    # use phylo to plot tree
    Phylo.draw(tree, 
               axes=axes[0], 
               do_show=False,
               label_colors=color_dict,
               )

    # remove artists
    axes[0].spines[['right', 'top', "left"]].set_visible(False)
    axes[0].set_yticks([])
    axes[0].set_ylabel("")

    # remove texts
    strain_pos = {}
    for txt in axes[0].texts:

        strain = txt.get_text().strip(" ")

        # remove the nodes
        if strain not in colors_all["strain"]:
            txt.set_visible(False)

        # get coordinates of texts
        else:
            txt.set_text(f" {strain} ")
            strain_pos[strain] = txt.get_position()[1]


    # get the y-lim
    ylim = axes[0].get_ylim()
    xlim = axes[0].get_xlim()
    axes[0].set_xlim(xlim[0], xlim[1]*1)
    axes[0].set_xlabel("Branch length")



    ### ADD THE POINT PLOT WITH THE RECOVERED STRAINS
    evolved_coms = evolved_sc["syncom"].unique(maintain_order=True)
    sc_ls = []
    for evolvedCom in evolved_coms:
        
        # subset
        evolvedComdf = evolved_sc.filter(pl.col("syncom") == evolvedCom).with_columns(pl.col("strain").alias("Evolved"))
        sc_ls.append((pl.DataFrame({"strain": strain_pos.keys()}).
                        with_columns(pl.lit(evolvedCom).alias("syncom"))
                        .join(colors_all.select("strain", "family"), how="left", on="strain")
                        .with_columns(pl.lit(evolvedComdf["condition"].unique()[0]).alias("condition"))
                        .join(evolvedComdf.select("strain", "Evolved"), how="left", on="strain")
                        .fill_null("No")))
    sc_df = pl.concat(sc_ls)

    # shaded background for non-native conditions
    axes[1].set_clip_on(False)
    for e in range(evolved_coms.len()):
        if (e+1)%2==1:
            axes[1].axvspan(-0.5, 0.5, facecolor='lightgrey', alpha=0.5)
            axes[1].text(x=e+0.5,
                         y=-1,
                         s=evolved_coms[e].split("-")[1],
                         ha="center"
                         )

    # plot the scatter
    color_dict["No"] = "k"
    sizes = {label: 5 if label == "No" else 50 for label in color_dict.keys()}
    sns.scatterplot(data=sc_df.to_pandas(), 
                            x="syncom", 
                            y="strain", 
                            ax=axes[1], 
                            hue="Evolved",
                            palette=color_dict,
                            legend=None,
                            size="Evolved",
                            sizes=sizes,
                            linewidth=0,
                            style="condition",
                            markers=markers,
                            )
    
    # move axis to right
    axes[1].yaxis.set_visible(False)
    axes[1].set_ylabel("")
    axes[1].spines[["left", "right"]].set_visible(False)
    axes[1].set_xlim(-0.5, 1.5)
    axes[1].set_xlabel("")

    # format xticklabels
    xticklabels = []
    for text in axes[1].get_xticklabels():
        text.set_text(title_dict[text.get_text().split("-")[0]])
        xticklabels.append(text)
    axes[1].set_xticklabels(xticklabels, ha="right", rotation=30, rotation_mode="anchor")



    ### MAKE GENOMIC BARPLOTS SHOWING LOCATION OF VARIANTS
    offset_x = pl.DataFrame({"strain": strain_pos.keys()}).with_columns(pl.lit(0).alias("offset_x"))
    for contig in genome_sizes["contig"].unique(maintain_order=True):

        # subset the data
        contig_subset = (genome_sizes
                         
                            # filtering of respective contigs
                            .filter(pl.col("strain").is_in(strain_pos.keys()))
                            .filter(pl.col("contig") == contig) 

                            # applying the offset on xaxis for stacked plots
                            .join(offset_x, how="left", on="strain") 
                            .with_columns(pl.col("offset_x").alias("left"))
                            .with_columns(pl.col("endpos").alias("width"))

                            # add y-pos and put into order
                            .join(pl.DataFrame({"strain": strain_pos.keys(), "y": strain_pos.values()}), how="left", on="strain")
                            .sort(pl.col("strain").cast(pl.Enum(list(strain_pos.keys()))))
                            )


        # subset variants
        variant_subset = (variants
                            .filter(pl.col("syncom") == syncom)
                            .filter(pl.col("contig") == contig)
                            .join(pl.DataFrame({"strain": strain_pos.keys(), "y": strain_pos.values()}), how="left", on="strain")

                            # applying the offset on xaxis for stacked plots
                            .join(offset_x, how="left", on="strain") 
                            .with_columns( (pl.col("loc") + pl.col("offset_x")).alias("loc") )
                            .filter(pl.col("y").is_not_null())
                            
                            # add color
                            .join(tree_colors, how="left", on=["syncom", "ID"])
                            .filter(pl.col("significant")) 
                            .filter(pl.col("color").is_not_null())
                            )
        
        # adjust the alpha
        variant_subset = variant_subset.select("loc", "y", "color", "edgecolor", "AF", "ID").unique()
        


        if contig_subset.shape[0]>0:
            axes[2].barh(y=contig_subset["y"].to_numpy(),
                         left=contig_subset["left"].to_numpy(),
                         width=contig_subset["width"].to_numpy(),
                         color="lightgrey",
                         edgecolor="lightgrey",
                        )
            
            if variant_subset.shape[0]>0:
                # create path collections
                af = variant_subset["AF"].to_numpy()
                face_rgba = np.array([mcolors.to_rgba(c, alpha=a) for c, a in zip(variant_subset["color"], af) if c != "lightgrey"])
                edge_rgba = np.array([mcolors.to_rgba(c, alpha=a) for c, a in zip(variant_subset["edgecolor"], af) if c == "black"])
                rect_faces = [Rectangle(xy=(row["loc"]-1e-2, row["y"]-0.425), width=4e-2, height=0.85) for row in variant_subset.iter_rows(named=True) if row["color"]!="lightgrey"]
                rect_edges = [Rectangle(xy=(row["loc"]-2e-2, row["y"]-0.41), width=4e-2, height=0.82) for row in variant_subset.iter_rows(named=True) if row["color"]=="lightgrey"]
                pc_faces = PatchCollection(rect_faces, facecolor=face_rgba, edgecolor="none")
                pc_edges = PatchCollection(rect_edges, facecolor="none", edgecolor=edge_rgba, linewidth=0.5)
                axes[2].add_collection(pc_faces)
                axes[2].add_collection(pc_edges)


            # now, update the offset_x locations
            offset_x = (contig_subset
                            .select("strain", "left", "width")
                            .with_columns((pl.col("left")+pl.col("width")+0.1).alias("offset_x"))
                            .drop("left", "width")
                       )

        
    # modify xticks
    axes[2].yaxis.set_visible(False)
    axes[2].spines[['right', 'top', "left"]].set_visible(False)
    axes[2].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}Mb"))

    # define the locators
    axes[2].xaxis.set_major_locator(ticker.FixedLocator([0,2,4,6,8,10]))

    # handle axis's
    axes[2].set_xlabel("Genomic location")
    axes[2].set_ylim(ylim)
    xlim_map = axes[2].get_xlim()
    axes[2].set_xlim(-0.05, xlim_map[1])


    # add legend here
    if syncom == "at":
        legend_color = [Patch(label="Evolved host", alpha=0)] + [Patch(facecolor=color if color!="black" else "white", edgecolor=color, label=specificity, alpha=1) for specificity, color in colors.items()]
        lgd = axes[2].legend(handles=legend_color, 
                        frameon=False,
                        )

        for item, label in zip(lgd.legend_handles, lgd.texts):
            if label._text  in ["Evolved host"]:
                width=item.get_window_extent(fig.canvas.get_renderer()).width
                label.set_position((-1.5*width,0))
                label.set_fontsize(12)



    # now manually adjust the subplot positions
    # x0, y0, width , heigh
    branch_width = 0.15
    syncom_width = 0.025*n_syncoms
    syncom_start = 0.125+branch_width + 0.05
    map_start=syncom_start + syncom_width+0.015
    for idx, ax in enumerate(fig.axes):
        
        # Bbox([[0.125, 0.10999999999999999], [0.33540723981900455, 0.88]])
        if idx == 0:
            ax.set_position((0.125, 0.1, branch_width, 0.88))

        # Bbox([[0.38099547511312215, 0.10999999999999999], [0.43359728506787326, 0.88]])
        if idx == 1:
            ax.set_position((syncom_start, 0.1, syncom_width, 0.88))
        
        # Bbox([[0.4791855203619909, 0.10999999999999999], [0.8999999999999999, 0.88]])
        if idx == 2:
            ax.set_position((map_start, 0.1, 1-branch_width-syncom_width-0.05, 0.88))


    # save tree
    fig.savefig(f"figures/shotgun_tree_condensed_{syncom}.pdf",  bbox_inches="tight")
    fig.savefig(f"figures/shotgun_tree_condensed_{syncom}.png", dpi=600,  bbox_inches="tight")
    fig.clf()




