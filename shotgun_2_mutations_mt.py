
# script to replicate figure 2b from 
# Lenksi LTEE shotgun paper
# https://www.nature.com/articles/nature24287/figures/2


# import modules
import polars as pl
import matplotlib.pyplot as plt
import seaborn.objects as so
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import ticker
import warnings

# ~~~~~~~~~~~~ CONFIG ~~~~~~~~~~~~ #
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "font.family": "Helvetica"
})


# ~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *


# ~~~~~~~~~~~~ MAPPING ~~~~~~~~~~~~ #
mapping = load_mapping()


# and import the inputs
variants_long = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv", 
                            separator="\t",
                            schema_overrides={"line": str}
                            )
                            .filter(pl.col("AF") > 0)
                            .with_columns(("g" + pl.col("generation").cast(str) + "_l" + pl.col("line")).alias("sampleID"))
                            .filter(pl.col("type")!="upstream")
                )


#1 get number of mutations per line per timepoint
df_cols = ["line", "generation", "treatment", "ID", "AF"]
group_cols = ["line", "generation", "treatment"]
mt_df = (variants_long
            .select(df_cols) 
            .group_by(group_cols)
            .agg(pl.col("AF").sum().alias("mt"))
            .sort("line", "generation")
            )

fixed = (variants_long
            .select(df_cols) 
            .filter(pl.col("AF") >= 0.95)
            .group_by(group_cols)
            .agg(pl.col("AF").sum().alias("fixed"))
            .sort("line", "generation")
            )

mt_df = (mt_df
            .join(fixed, how="left", on=["line", "generation", "treatment"])
            .fill_null(0)
            )

# now the plotting
# x-axis will be time, y-axis will be mutations
# coloring will be according to strain | replicate line

# define colors
colors = {
    "1": "darkred",
    "2": "red",
    "3": "orangered",
    "4": "coral",
    "5": "salmon",
    "6.2": "black",
    "7.2": "dimgrey",
    "8.2": "slategrey",
    "9": "darkgrey",
    "10": "lightgrey",
    "11": "purple",
    "12": "darkviolet",
    "13": "blueviolet",
    "14": "magenta",
    "15": "violet",
    "16": "navy",
    "17": "blue",
    "18": "cornflowerblue",
    "19": "lightskyblue",
    "20": "deepskyblue",
}

names = {
    "1": "$At$-SC$^{At}$-1",
    "2": "$At$-SC$^{At}$-2",
    "3": "$At$-SC$^{At}$-3",
    "4": "$At$-SC$^{At}$-4",
    "5": "$At$-SC$^{At}$-5",
    "6.2": "$Lj$-SC$^{At}$-6",
    "7.2": "$Lj$-SC$^{At}$-7",
    "8.2": "$Lj$-SC$^{At}$-8",
    "9": "$Lj$-SC$^{At}$-9",
    "10": "$Lj$-SC$^{At}$-10",
    "11": "$At$-SC$^{Lj}$-11",
    "12": "$At$-SC$^{Lj}$-12",
    "13": "$At$-SC$^{Lj}$-13",
    "14": "$At$-SC$^{Lj}$-14",
    "15": "$At$-SC$^{Lj}$-15",
    "16": "$Lj$-SC$^{Lj}$-16",
    "17": "$Lj$-SC$^{Lj}$-17",
    "18": "$Lj$-SC$^{Lj}$-18",
    "19": "$Lj$-SC$^{Lj}$-19",
    "20": "$Lj$-SC$^{Lj}$-20",
}

markers = {"native": "s", "non-native": "D"}
linestyles = {"native": "-", "non-native": ":"}


fig, axes = plt.subplots(figsize=(5,3))
for treatment in mt_df["treatment"].unique():
    (so.Plot(data=mt_df.filter(pl.col("treatment") == treatment).to_pandas(),
            x="generation",
            y="mt",
            color="line",
            marker="treatment",)
            
            # add the artists
            .add(so.Line(linestyle=linestyles[treatment], pointsize=5, edgewidth=0, alpha=0.5))

            # scale artists
            .scale(color=colors)
            .scale(marker=markers)

            .on(axes)
            ).plot()


# add labels
axes.set_ylabel("Mutations, $M(t)$")
axes.set_xlabel("Plant cycle, $t$")

# reformat axis
axes.set_xticks(range(1,17))


# manage legend
[legend.set_visible(False) for legend in fig.legends]

# write
axes.spines[["top", "right"]].set_visible(False)
fig.tight_layout(pad=0)
fig.savefig(f"figures/shotgun_mutation_mt_line.png", dpi=600)
fig.savefig(f"figures/shotgun_mutation_mt_line.pdf")



# now vs fixed
fig, axes = plt.subplots(figsize=(5,3))
for treatment in mt_df["treatment"].unique():
    (so.Plot(data=mt_df.filter(pl.col("treatment") == treatment).to_pandas(),
            x="mt",
            y="fixed",
            color="line",
            marker="treatment",)
            
            # add the artists
            .add(so.Line(linestyle=linestyles[treatment], 
                         alpha=0.5, 
                         edgewidth=0, 
                         pointsize=0))

            # scale artists
            .scale(color=colors)
            .scale(marker=markers)

            .on(axes)
            ).plot()
    
    for generation in mt_df["generation"].unique():
        (so.Plot(data=mt_df.filter(pl.col("treatment") == treatment).filter(pl.col("generation") == generation).to_pandas(),
                 x="mt",
                 y="fixed",
                 color="line",
                 marker="treatment",
                )
                
                # add the artists
                .add(so.Dot(alpha=0.5, 
                            edgewidth=0,
                            pointsize=(generation**(1/1.5))))

                # scale artists
                .scale(color=colors)
                .scale(marker=markers)

                .on(axes)
                ).plot()



# add labels
axes.set_ylabel("Fixed mutations (AF ≥ 0.95)")
axes.set_xlabel("Mutations, $M(t)$")
axes.yaxis.set_major_locator(ticker.MaxNLocator(nbins=10, steps=[1,2,5,10,], integer=True, min_n_ticks=1))

# manage legend
[legend.set_visible(False) for legend in fig.legends]
axes.axline((0,0), (1,1), color="k", linestyle="--", alpha=0.5)

# write
axes.spines[["top", "right"]].set_visible(False)
fig.tight_layout(pad=0)
fig.savefig(f"figures/shotgun_mutation_mt_vs_fixed_line.png", dpi=600)
fig.savefig(f"figures/shotgun_mutation_mt_vs_fixed_line.pdf")







# make dedicated legend
fig, axes = plt.subplots(figsize=(3.2,3.2))
axes.xaxis.set_visible(False)
axes.yaxis.set_visible(False)
axes.spines[["top", "bottom", "left", "right"]].set_visible(False)

legend_color = [Patch(label="Population", alpha=0)] + [Patch(facecolor=colors[line], label=names[line], alpha=1) for line in colors.keys()]
legend_marker = [Patch(label="Condition", alpha=0)] + [Line2D([0], [0], marker=markers[label], linestyle=linestyles[label], label=label.title(), markerfacecolor='k', color="k") for label in markers.keys()]


lgd = axes.legend(handles=legend_color + legend_marker, 
                  loc='center', 
                  bbox_to_anchor=(0.5, 0.5),
                  ncol=2,
                  frameon=False
                  )

for item, label in zip(lgd.legend_handles, lgd.texts):
    if label._text  in ["Population", "Plant", "SynCom", "Generation", "Condition", "Enriched host"]:
        width=item.get_window_extent(fig.canvas.get_renderer()).width
        label.set_position((-1.5*width,0))


fig.savefig("figures/shotgun_mutations_mt_legend.pdf")
fig.savefig("figures/shotgun_mutations_mt_legend.png")


