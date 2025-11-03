
# goal is to illustrate the non-randomness of mutations, preferenetially using upset plots 
# or venn diagrams

# import modules
import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import warnings
import warnings
from matplotlib.patches import Wedge
from mpl_toolkits.axes_grid1.inset_locator import inset_axes




# ~~~~~~~~~~~~ CONFIGS ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "font.family": "Helvetica"
})

# ~~~~~~~~~~~~ MAPPING ~~~~~~~~~~~~ #
mapping = load_mapping()
kegg = load_kegg().drop("desc").group_by("ko").agg(pl.col(col).str.join(", ").alias(col) for col in ["c", "b", "a"])


# ~~~~~~~~~~~~ MAPPING ~~~~~~~~~~~~ #
level = "geneID"


# and import the inputs
variants_long = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv", 
                             separator="\t",
                             schema_overrides={"line": str}
                            )
                            .with_columns((pl.col("syncom") + "-" + pl.col("plant") + "-" + pl.col("line")).alias("line_name"))
                            .filter(pl.col(level).is_not_null())
                            .with_columns(pl.col("line").str.replace(".2", "", literal=True))
                )



def get_title(text):
    sc, plant, replicate = text.split("-")
    plant = "Lj" if plant == "gifu" else "At"
    return f"${sc.title()}$-SC-{replicate}$^" + "{" + plant + "}$"

def get_color(text):
    sc, plant, replicate = text.split("-")
    return "red" if plant == "col0" else "blue"


def draw_pie(ax, ratios, X, Y, size):
    # Create a small inset axes at (X,Y) in data coords
    inset = inset_axes(ax, width=size, height=size, 
                       loc='center',
                       bbox_to_anchor=(X, Y),
                       bbox_transform=ax.transData,
                       borderpad=0)
    start = 0
    for frac, color in ratios:
        end = start + frac * 360
        inset.add_patch(Wedge((0.5, 0.5), 0.5, start, end,
                              facecolor=color, edgecolor="white"))
        start = end
    inset.set_xlim(0,1)
    inset.set_ylim(0,1)
    inset.axis("off")



color_df_ls = []
gene_counter = {"single-hit": 0, "multi-hit": 0}
system_adaptation = []
for syncom in ["at", "lj"]:



    # get ids and sunset data
    variant_subset = (variants_long
                        .filter(pl.col("syncom") == syncom)
                        .with_columns((pl.col("syncom") + "-" + pl.col("plant") + "-" + pl.col("line")).alias("line_name")))
    

    # get the total counts
    counts_line = {}
    for line in variant_subset.sort("syncom", "plant", pl.col("line").cast(float), descending=True)["line_name"].unique(maintain_order=True):
        replicate_subset = variant_subset.filter(pl.col("line_name") == line)
        counts_line[f"{line}"] = len(replicate_subset[level].unique().to_list())


    # now get the insertion sizes, but the ones that are aggregated into groups
    insertions = {}
    for level_val in variant_subset[level].unique():
        level_subset = variant_subset.filter(pl.col(level) == level_val)
        insertions[level_val] = level_subset["line_name"].unique().sort().to_list()
    

    # get the solid system adaptation
    system_adaptation += [key for key, value in insertions.items() if len(value)==10]

    # we have 4 color categories:
    # grey - line specific
    # red - arabidopsis specific
    # blue - lotus specific
    # black - unspecific
    # now populate these:
    color_d = {}
    for level_val in insertions.keys():
        if len(insertions[level_val]) == 1:
            color_id = "grey"
            gene_counter["single-hit"] += 1
        elif all(["col0" in e for e in insertions[level_val]]):
            color_id = "red"
            gene_counter["multi-hit"] += 1
        elif all(["gifu" in e for e in insertions[level_val]]):
            color_id = "blue"
            gene_counter["multi-hit"] += 1
        else:
            color_id = "black"
            gene_counter["multi-hit"] += 1
        color_d[level_val] = color_id

    # here save a table with the color mappings, containing at least the following columns
    # ID ┆ syncom ┆ enriched (plant)
    color_df = (pl.DataFrame({"syncom": syncom, "geneID": color_d.keys(), "color": color_d.values()})
                    .join(variants_long.select("syncom", "geneID", "ID").unique(), how="left", on=["syncom", "geneID"])
                    .sort("geneID")
               )
    color_df_ls.append(color_df)

    
    # now generate the dots of the upset plot
    groups_d = {member: 0 for member in variant_subset["line_name"].unique()}
    groups_d["col0"] = 0
    groups_d["gifu"] = 0
    groups_d["unspecific"] = 0
    for member in insertions.values():
        if len(member) == 1:
            groups_d[member[0]] += 1
        elif all(["col0" in e for e in member]):
            groups_d["col0"] += 1
        elif all(["gifu" in e for e in member]):
            groups_d["gifu"] += 1
        else:
            groups_d["unspecific"] += 1
    groups_d = dict(sorted(groups_d.items(), key=lambda item: item[1], reverse=True))

    # determine colors
    color_mask = np.array(list(groups_d.keys()))
    color_vec = np.full(color_mask.shape, "grey", dtype="<U20")
    color_vec[color_mask == "col0"] = "red"
    color_vec[color_mask == "gifu"] = "blue"
    color_vec[color_mask == "unspecific"] = "black"


    # start plotting here
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8,5),height_ratios=(1,1.5))


    # plot the intersections
    axes[0].bar(x=groups_d.keys(), height=groups_d.values(), color=color_vec)
    axes[0].xaxis.set_visible(False)
    axes[0].margins(x=0.025, y=0)
    axes[0].spines[["right", "top"]].set_visible(False)
    axes[0].set_ylim(0, (np.max(list(groups_d.values()))*1.1))
    axes[0].set_ylabel("Intersection size")
    axes[0].set_title(f"Mutational overlap in ${syncom.title()}$-SC populations (ORFs)")

    # texts
    for key, val in groups_d.items():
        axes[0].text(x=key, 
                       y=val, 
                       s=val,
                       ha="center",
                       va="bottom",
                       color="k")

    # now populate the dot matrix
    axes[1].scatter(x=list(groups_d.keys())*len(counts_line.keys()), y=list(counts_line.keys())*len(groups_d.keys()), color="lightgrey")
    axes[1].spines[["top", "bottom", "right"]].set_visible(False)
    line_d = {}
    for group_idx, group in enumerate(groups_d.keys()):
        if group in ["col0", "gifu", "unspecific"]:
            line_d[group] = []
        for line_idx, line in enumerate(counts_line.keys()):

            # plot solid grey points for specific variants
            if line == group:
                draw_pie(axes[1], [(1, "grey")], group_idx, line_idx, 0.15)

            # now plot the arabidopsis points
            elif ("col0" in line) and (group == "col0"):
                genesIDs = color_df.filter(pl.col("color")=="red")["geneID"].unique()
                total = len(genesIDs)
                in_line = 0
                for geneID in genesIDs:
                    if line in insertions[geneID]:
                        in_line+=1
                fraction = in_line/total
                draw_pie(axes[1], [(fraction, "red")], group_idx, line_idx, 0.25)
                line_d[group].append(line_idx)

            elif ("gifu" in line) and (group == "gifu"):
                genesIDs = color_df.filter(pl.col("color")=="blue")["geneID"].unique()
                total = len(genesIDs)
                in_line = 0
                for geneID in genesIDs:
                    if line in insertions[geneID]:
                        in_line+=1
                fraction = in_line/total
                draw_pie(axes[1], [(fraction, "blue")], group_idx, line_idx, 0.25)
                line_d[group].append(line_idx)

            elif group == "unspecific":
                genesIDs = color_df.filter(pl.col("color")=="black")["geneID"].unique()
                total = len(genesIDs)
                in_line = 0
                for geneID in genesIDs:
                    if line in insertions[geneID]:
                        in_line+=1
                fraction = in_line/total
                draw_pie(axes[1], [(fraction, "black")], group_idx, line_idx, 0.25)
                line_d[group].append(line_idx)
    
    shorting = 0.35
    for group_idx, group in enumerate(line_d.keys()):
        y_vals = [[a+shorting, b-shorting] for a, b in zip(line_d[group], line_d[group][1:])]
        for y in y_vals:
            axes[1].plot([group]*len(y), y, color="red" if group == "col0" else "blue" if group == "gifu" else "black")
    axes[1].margins(x=0.06, y=0.05)
    axes[1].xaxis.set_visible(False)


    # reformat yticks
    tick_list = []
    for xtick in axes[1].get_yticklabels():
        condi = xtick.get_text()
        xtick.set_color(get_color(condi))
        xtick.set_text(get_title(condi))
        tick_list.append(xtick)
    axes[1].set_yticklabels(tick_list)

    # save figure
    fig.tight_layout(h_pad=0.1)
    fig.savefig(f"figures/shotgun_upset_plot_{syncom}_custom.pdf")
    fig.savefig(f"figures/shotgun_upset_plot_{syncom}_custom.png")


# aggregate the table with colors
color_df = pl.concat(color_df_ls)
color_df.write_csv("results/tree_colors.tsv", separator="\t")


# print multi-hit genes
print(f"number of multi-hit genes: {gene_counter["multi-hit"]} vs. single-hit {gene_counter["single-hit"]}")


# write table with system adaptation
system_adaptation = (pl.DataFrame({"geneID": system_adaptation}).join(loadGff()
                                                                      .select("geneID", "ko", "alias", "desc").unique(), 
                                                                      how="left", 
                                                                      on="geneID"))
system_adaptation.write_csv("results/system_adaptation.tsv", separator="\t")




# just generate the legend
colors = {"$A. thaliana$ (Col-0)": "red",
          "$L. japonicus$ (Gifu)": "blue",
          "Trajectory": "grey",
          "Unspecific": "black"
          }
legend_color = [Patch(label="Mutation specificity", alpha=0)] + [Patch(facecolor=color, label=specific, alpha=1) for specific, color in colors.items()]

fig, axes = plt.subplots(figsize=(2.5, 1), 
                        )

lgd = axes.legend(handles=legend_color, 
                  loc='center', 
                  bbox_to_anchor=(0.5, 0.5),
                  frameon=False,
                  )
axes.xaxis.set_visible(False)
axes.yaxis.set_visible(False)
axes.spines[["top", "bottom", "left", "right"]].set_visible(False)

for item, label in zip(lgd.legend_handles, lgd.texts):
    if label._text  in ["Mutation specificity", "SynCom", "Generation", "Condition", "Enriched host"]:
        width=item.get_window_extent(fig.canvas.get_renderer()).width
        label.set_position((-1.5*width,0))
        label.set_fontsize(12)

fig.tight_layout()
fig.savefig(f"figures/shotgun_upset_plot_legend.png", dpi=600)
fig.savefig(f"figures/shotgun_upset_plot_legend.pdf")
