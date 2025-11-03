
# script to plot the volcano hits


# import modules
import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import seaborn.objects as so
import textalloc as ta
import warnings
warnings.filterwarnings("ignore")




# ~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *
plt.rcParams.update({
    "font.family": "Helvetica"
})

# ~~~~~~~~~~~~ MAPPING ~~~~~~~~~~~~ #
mapping = load_mapping()
kegg_tab = load_kegg()
gff = loadGff()


# ~~~~~~~~~~~~ TABLE ~~~~~~~~~~~~ #
glm_result_df_meta = pl.read_csv("results/volcano.tsv", separator="\t").with_columns(pl.col("enriched").str.replace("NS", "n.s."))


# figure aesthetics
colors = {"col0": "red", "gifu": "blue", "n.s.": "grey"}
markers = {"native": "s", "non-native": "D", "n.s.": "X"}



# plot two vulcano plots
for idx, syncom in enumerate(["at", "lj"]):

    # make 2 figures
    fig, axes = plt.subplots(figsize=(5,5))

    #subset the syncom data to plot
    sc_subset = glm_result_df_meta.filter(pl.col("syncom") == syncom)

    # get number of significant genes
    n_significant = sc_subset.filter(pl.col("significant")).group_by("enriched").len().sort("enriched")

    (so.Plot(data = sc_subset.to_pandas(),
             x="ratio_plant",
             y="-log10(pval)",
             color="enriched",
             alpha="annotated",
             )
             
             # add the artists
             .add(so.Dot(marker="o", edgewidth=0))

             # define axes
             .on(axes)

             # add colors manualy
             .scale(color=colors)
             .scale(alpha={True: 1, False: 0.5})

             #draw
             .plot()
             )
    
    # add the gene annotations
    for enriched in ["col0", "gifu"]:
        sc_subset_ = (sc_subset
                        .filter(pl.col("syncom") == syncom)
                        .filter(pl.col("enriched")==enriched)
                        .select(["syncom", "enriched", "geneID", "alias"])
                        .join(sc_subset.select(["syncom", "enriched", "geneID", "ratio_plant", "-log10(pval)"]).unique(), how="left", on=["syncom", "enriched", "geneID"])
                     )
        try:
            sc_subset_ = pl.concat(
                df.sample(n=1) for df in 
                sc_subset_.partition_by(["geneID"], include_key=True)
            )
            texts = sc_subset_["alias"].to_numpy()
            ta.allocate(ax = axes,
                        x=sc_subset_["ratio_plant"].to_numpy()[texts != None],
                        y=sc_subset_["-log10(pval)"].to_numpy()[texts != None],
                        text_list=texts[texts != None],
                        x_scatter=sc_subset_["ratio_plant"].to_numpy()[texts != None],
                        y_scatter=sc_subset_["-log10(pval)"].to_numpy()[texts != None],
                        textsize=8,
                        nbr_candidates=30,
                        linecolor="k",
                        draw_all=False,
                        min_distance=0,
                        direction="north",
                        **{"ha": "center"}
                        )
        except:
            pass

    axes.set_title(f"Significant variants for ${syncom.title()}$-SC")
    axes.set_yscale('symlog')
    axes.set_ylim((0, sc_subset["-log10(pval)"].max()*2))
    xlims = axes.get_xlim()
    offset = np.abs(xlims).sum()*0.1
    axes.set_xlim((xlims[0]-offset, xlims[1]+offset))
    axes.set_xlabel("Log2 FC  $A. thaliana$ vs $L. japonicus$")
    axes.spines[["top", "right"]].set_visible(False)

    # add dashed lines indicating the cutoff for the significance
    axes.axvline(x=-1, color="k", linestyle=":", alpha=0.5)
    axes.axvline(x=1, color="k", linestyle=":", alpha=0.5)
    axes.axhline(y=-np.log10(0.05), color="k", linestyle=":", alpha=0.5)

    # add text with number of significant hits
    x_axis_range = np.sum(np.abs(axes.get_xlim()))*0.025
    x_offset_text = [x_axis_range, -x_axis_range]
    for idx_plant, row in enumerate(n_significant.iter_rows()):
        axes.text(x=axes.get_xlim()[idx_plant]+x_offset_text[idx_plant],
                  y=-np.log10(0.05)*1.1,
                  s=f"$n={row[1]}$",
                  va="center",
                  ha="left" if idx_plant == 0 else "right")
    
    condition = ["native", "non-native"] if syncom == "at" else ["non-native", "native"]
    for idx_plant, row in enumerate(n_significant.iter_rows()):
        axes.text(x=(axes.get_xlim()[idx_plant]+(1 if idx_plant==1 else -1))/2,
                  y=axes.get_ylim()[1]*0.85,
                  s=condition[idx_plant].title(),
                  va="center",
                  ha="center")


    # manage legend
    [legend.set_visible(False) for legend in fig.legends]

    # some figure aestetics
    axes.set_ylabel("-$Log$$_{10}$($P$-value)")
    fig.tight_layout()


    # write
    fig.savefig(f"figures/shotgun_vulcano_all_ID_{syncom}.png", dpi=600)
    fig.savefig(f"figures/shotgun_vulcano_all_ID_{syncom}.pdf")



# just generate the legend
legend_color = [Patch(label="Enriched host", alpha=0)] + [Patch(facecolor=color, label=get_label_enriched(enriched), alpha=1) for enriched, color in colors.items()]

fig, axes = plt.subplots(figsize=(2, 1), 
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
    if label._text  in ["Plant", "SynCom", "Generation", "Condition", "Enriched host"]:
        width=item.get_window_extent(fig.canvas.get_renderer()).width
        label.set_position((-1.5*width,0))
        label.set_fontsize(12)

fig.tight_layout()
fig.savefig(f"figures/shotgun_vulcano_all_ID_legend.png", dpi=600)
fig.savefig(f"figures/shotgun_vulcano_all_ID_legend.pdf")

