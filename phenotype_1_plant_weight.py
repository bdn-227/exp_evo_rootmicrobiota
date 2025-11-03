
# script by niklas kiel to analyse weight of plants 
# during the LTEE

import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import false_discovery_control, mannwhitneyu, spearmanr

plt.rcParams.update({
    "font.family": "Helvetica"
})



# load the plant weight
weight = pl.read_csv("raw_data/plant_weights.tsv", separator="\t", schema_overrides={"line": str})


# first, plot weights for evolution
weight_axenic = weight.filter(pl.col("sample_type") == "axenic").select("generation", "plant", "weight").group_by("generation", "plant").agg(pl.col("weight").mean().alias("axenic"))
weight_evo = weight.filter(pl.col("sample_type") == "generation").with_columns((pl.col("generation")+1)).join(weight_axenic, how="left", on=["plant", "generation"]).with_columns((pl.col("weight")/pl.col("axenic")).alias("weight"))


def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "NS"

def number_formatter(number, precision=3):
    """
    function to reformat numbers into nice scientific notation
    """
    
    if "e" in str(number):
        return np.format_float_scientific(number, precision=precision).replace("e", "×$10^{") + "}$"
    else:
        return np.round(number, precision)

# create legend
colors = {"at": "#eca6b5", "lj": "#90d8f5"}
markers = {"Native": "s", "Non-native": "D"}
legend_color = [Patch(label="SynCom", alpha=0)] + [Patch(facecolor=color, label=f"${label.title()}$-SC", alpha=1) for label, color in colors.items()]
legend_marker = [Patch(label="Condition", alpha=0)] + [Line2D([0], [0], marker=marker, linestyle=None, label=label, markerfacecolor='k', color="w") for label, marker in markers.items()]



for plant in ["col0", "gifu"]:
    plant_subset = (weight_evo
                        .filter(pl.col("plant") == plant)
                        .with_columns(pl.when(pl.col("syncom") == "at").then(pl.col("generation")-0.15)
                                        .when(pl.col("syncom") == "lj").then(pl.col("generation")+0.15)
                                        .otherwise("generation").alias("x"))
                        .with_columns(pl.when(pl.col("syncom") == "at").then(pl.lit("#eca6b5"))
                                        .when(pl.col("syncom") == "lj").then(pl.lit("#90d8f5"))
                                        .otherwise(pl.lit("grey")).alias("color"))
                        .with_columns(pl.when(pl.col("condition") == "native").then(pl.lit("s"))
                                        .when(pl.col("condition") == "non-native").then(pl.lit("D"))
                                        .otherwise(pl.lit("grey")).alias("marker"))
                        
                        )
    
    
    # calculate statistics
    sc_ls = ["at", "lj"]
    corr_sc_ls=[]
    corr_pval_ls = []
    corr_ls = []

    comp_sc_ls=[]
    comp_pval_ls = []
    comp_x_ls = []
    comp_y_ls = []
    for syncom in sc_ls:
        sc_subset = plant_subset.filter(pl.col("syncom") == syncom)
        corr_res = spearmanr(sc_subset["generation"].to_numpy(), 
                             sc_subset["weight"].to_numpy())
        corr_pval_ls.append(corr_res.pvalue)
        corr_ls.append(corr_res.statistic)
        corr_sc_ls.append(syncom)


        initial = sc_subset.filter(pl.col("generation") == 1)
        for generation in np.arange(sc_subset["generation"].min()+1, sc_subset["generation"].max()+1):
            gen_subset = sc_subset.filter(pl.col("generation") == generation)
            comp_res = mannwhitneyu(x=initial["weight"],
                                    y=gen_subset["weight"])
            comp_sc_ls.append(syncom)
            comp_pval_ls.append(comp_res.pvalue)
            comp_x_ls.append(gen_subset["x"].mean())
            comp_y_ls.append((gen_subset["weight"].mean() + gen_subset["weight"].std())*1.1)
        
    corr_padj_ls = false_discovery_control(corr_pval_ls)
    comp_padj_ls = [convert_pvalue_to_asterisks(padj) for padj in false_discovery_control(comp_pval_ls)]


    # plot that shit
    fig, axes = plt.subplots(figsize=(10,3))

    # make point plot with means
    print(f"plant: {plant} - n = {plant_subset.shape[0]}")
    for marker in plant_subset["marker"].unique():
        axes.scatter(x=plant_subset.filter(pl.col("marker")==marker)["x"].to_numpy(),
                     y=plant_subset.filter(pl.col("marker")==marker)["weight"].to_numpy(),
                     color=plant_subset.filter(pl.col("marker")==marker)["color"].to_numpy(),
                     marker=marker,
                     alpha=0.1,
                     linewidths=0,
                     zorder=5
                    )

    # plot error bars
    stats = plant_subset.group_by("x", "color", "marker").agg(pl.col("weight").std().alias("std"), pl.col("weight").mean().alias("mean"))
    axes.errorbar(x=stats["x"].to_numpy(),
                  y=stats["mean"].to_numpy(),
                  yerr=stats["std"].to_numpy(),
                  ecolor=stats["color"].to_numpy(),
                  fmt = 'none',
                  zorder=4
                 )

    # plot means
    for marker in plant_subset["marker"].unique():
        axes.scatter(x=stats.filter(pl.col("marker")==marker)["x"].to_numpy(),
                     y=stats.filter(pl.col("marker")==marker)["mean"].to_numpy(),
                     color=stats.filter(pl.col("marker")==marker)["color"].to_numpy(),
                     marker=marker,
                     linewidths=0,
                     zorder=4
                    )
                    

    axes.set_ylim((0, 6))
        
    # write mann-withney results
    for sc, padj, x, y in zip(comp_sc_ls, comp_padj_ls, comp_x_ls, comp_y_ls):
        axes.text(x=x,
                  y=y,
                  s=padj,
                  ha="center",
                  color=colors[sc],
                  fontsize=8,
                  )


    # manage axes
    axes.set_xticks(np.arange(1,17))
    axes.set_xlim((0.5, 16.5))

    # labels
    axes.set_ylabel("Normalized weight")
    axes.set_xlabel("Plant cycle, $t$")
    axes.set_title("Plant weight of " + ("$A. thaliana$ (Col-0)" if plant == "col0" else "$L. japonicus$ (Gifu)") + " plants")

    # remove spines
    axes.spines[["top", "right"]].set_visible(False)

    # manage legend
    [legend.set_visible(False) for legend in fig.legends]
    lgd = axes.legend(handles=legend_color + legend_marker, 
                      loc='upper right', 
                      frameon=False,
                      ncol=2,
                      #bbox_to_anchor=(0.5, 0.5),
                     )
    for item, label in zip(lgd.legend_handles, lgd.texts):
        if label._text  in ["Plant", "SynCom", "Generation", "Condition"]:
            width=item.get_window_extent(fig.canvas.get_renderer()).width
            label.set_position((-1.5*width,0))

    # save figure
    fig.tight_layout(pad=0)
    fig.savefig(f"figures/weight_generation_{plant}.pdf")
    fig.savefig(f"figures/weight_generation_{plant}.png")
    fig.savefig(f"figures/weight_generation_{plant}.svg")





