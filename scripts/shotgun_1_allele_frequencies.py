
# ~~~~~~~~~~~~ IMPORTS ~~~~~~~~~~~~ #
import polars as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import seaborn.objects as so
from matplotlib.patches import Patch
from matplotlib import ticker
import warnings


# configs
mpl.use("Agg")
warnings.filterwarnings("ignore")
pl.Config.set_fmt_str_lengths(200)
pl.Config(tbl_rows=100)
plt.rcParams.update({
    "font.family": "Helvetica"
})

# custom functions
from shotgun_analysis_functions import *

# parameters
separate_by = "plant"
afcutoff=0


# ~~~~~~~~~~~~ LOAD DATA ~~~~~~~~~~~~ #
mapping = load_mapping()
gff = loadGff()


# load the variants
variants_long = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv",
                             separator="\t",
                             schema_overrides={"line": str}
                            )
                            .filter(pl.col("type")!="upstream")
                            .filter(pl.col("generation") > 0)
                            .sample(fraction=1.0, shuffle=True)
                            .with_columns((pl.col("ID_line") + "-" + pl.col("strain")).alias("id_strain"))
                            .join(load_colors(), how="left", on="strain")
                            )

# get index for the fixed ones
lines2highlight = variants_long.filter(pl.col("AF") >= 0.95)["id_strain"]

# get some stats for writing
variants_long["ID_line"].unique().len()
variants_long["strain"].unique().len()

variants_long.filter(pl.col("AF") >= 0.95)["ID"].unique().len()
variants_long.filter(pl.col("AF") >= 0.95)["strain"].unique().len()


# make a dict for line-stypes
markers = {"native": "s", "non-native": "D"}
linestyles = {"native": "-", "non-native": "--"}
colors = {id_strain: color for id_strain, color in zip(variants_long["id_strain"], variants_long["color"])}
title_dict = {
            'col0': "$A. thaliana$ (Col-0)", 
            'gifu': "$L. japonicus$ (Gifu)", 
            }

figsize=(10,2)
for syncom in ["at", "lj"]:

    # subset the variants table
    variants_wide_sc = variants_long.filter(pl.col("syncom") == syncom)
    print(f"total mutations found in {syncom}: " + str(variants_wide_sc["ID_line"].unique().len()))

    # create figure
    fig, axes = plt.subplots(nrows=1,
                             ncols=1, 
                             figsize=figsize,
                             constrained_layout=True,
                            )


    ############### ALLELE FREQUENCY PLOT ###############
    (so.Plot(data=variants_wide_sc.to_pandas(),
             x="generation",
             y="AF",
             color="id_strain",
             )
             
             .add(so.Line(alpha=0.1, linewidth=0.5, pointsize=0, color="grey"), legend=False)
             .on(axes)
             .scale(marker=markers)
             .scale(color=colors)

             ).plot()
    
    # plot padded values
    (so.Plot(data=variants_wide_sc.filter((pl.col("id_strain").is_in(lines2highlight))).filter(pl.col("af_type").is_in(["zero", "original"])).to_pandas(),
             x="generation",
             y="AF",
             color="id_strain",
             marker="treatment")
             
             .add(so.Line(alpha=0.5, linewidth=1.5, pointsize=1.5, edgewidth=0), legend=False)
             .on(axes)
             .scale(marker=markers)
             .scale(color=colors)

             ).plot()

    x_min = 1
    x_max = 16
    x_offset = 0.25

    # figure aesthetics
    y = 1.05
    axes.set_ylim([0, y])
    axes.set_xlim(1-x_offset, 16+x_offset)
    axes.set_xticks(range(1,17))

    # axis labels
    axes.set_ylabel("Allele frequency, $f(t)$")
    axes.set_xlabel("Plant cycle, $t$")
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    ############### ALLELE FREQUENCY PLOT ###############

    # save
    fig.suptitle(f"Allele frequencies of ${syncom.title()}$-SC", y=0.975)
    fig.tight_layout(pad=0)
    fig.savefig(f"figures/shotgun_allele_frequencies_{syncom}.pdf")
    fig.savefig(f"figures/shotgun_allele_frequencies_{syncom}.png", dpi=600)


# make a figure containing the legend for the community profiles over time
colors = color_dict_family()

# make custom legend
legend_color = [Patch(label="Family", alpha=0)] + [Patch(facecolor=color, label=label, alpha=0.6) for label, color in colors.items()]


# just generate the legend
fig, axs = plt.subplots(nrows=1,
                        ncols=1, 
                        figsize=(2, 3.8), 
                        )

lgd = axs.legend(handles=legend_color, 
                 loc='center', 
                 bbox_to_anchor=(0.5, 0.5),
                 ncol=1,
                 frameon=False,
                 )

axs.xaxis.set_visible(False)
axs.yaxis.set_visible(False)
axs.spines[["top", "bottom", "left", "right"]].set_visible(False)

for item, label in zip(lgd.legend_handles, lgd.texts):
    if label._text  in ["Family"]:
        label.set_position((-30,0))
        label.set_fontsize(12)

fig.tight_layout(pad=0)
fig.savefig(f"figures/shotgun_family_legend_vertical.png", dpi=600)
fig.savefig(f"figures/shotgun_family_legend_vertical.pdf")










############### DENSITY HISTOGRAM ###############
markers = {"native": "s", "non-native": "D"}
linestyles = {"native": "-", "non-native": "--"}
colors = {"col0": "red", "gifu": "blue"}




for syncom in ["at", "lj"]:
    
    # subset the data
    variants_wide_sc = variants_long.filter(pl.col("syncom") == syncom)

    # plotting
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5,3))
    for plant in ["col0", "gifu"]:
        axes.hist(variants_wide_sc.filter(pl.col("AF")>0).filter(pl.col("plant")==plant)["AF"].to_numpy(),
                  bins=25,
                  range=(0,1),
                  histtype='step',
                  orientation='vertical',
                  linestyle=linestyles[variants_wide_sc.filter(pl.col("plant")==plant)["treatment"].max()],
                  color=colors[plant],
                  linewidth=2,
                 )

    axes.set_xlabel("Allele frequency, $f(t)$")
    axes.set_ylabel("Count")
    axes.xaxis.set_major_locator(ticker.MultipleLocator(0.1))

    # remove borders
    axes.spines[['top', 'right']].set_visible(False)

    # titles and stuff
    axes.set_title(f"Allele frequency distributions of ${syncom.title()}$-SC", ha="center")
    axes.margins(x=0)
    axes.set_xlim((0,1.01))

    fig.tight_layout(pad=0)
    fig.savefig(f"figures/shotgun_allele_frequencies_hist_{syncom}.pdf")
    fig.savefig(f"figures/shotgun_allele_frequencies_hist_{syncom}.png", dpi=600)




# create legend for histogram
legend_color = [Patch(label="Plant", alpha=0)] + [Patch(facecolor=color, label=title_dict[label], alpha=1) for label, color in colors.items()]
legend_marker = [Patch(label="Condition", alpha=0)] + [Line2D([0], [0], linestyle=linestyles[label], label=label.title(), markerfacecolor='k', color="k") for label in markers.keys()]


fig, axes = plt.subplots(figsize=(2, 1.5))

lgd = axes.legend(handles=legend_color + legend_marker, 
                  loc='center', 
                  bbox_to_anchor=(0.5, 0.5),
                  frameon=False,
                  )
axes.xaxis.set_visible(False)
axes.yaxis.set_visible(False)
axes.spines[["top", "bottom", "left", "right"]].set_visible(False)

for item, label in zip(lgd.legend_handles, lgd.texts):
    if label._text  in ["Fitness on"]:
        width=item.get_window_extent(fig.canvas.get_renderer()).width
        label.set_position((-1.5*width,0))
        label.set_fontsize(12)

fig.tight_layout(pad=0)
fig.savefig(f"figures/shotgun_histogram_legend.pdf")
fig.savefig(f"figures/shotgun_histogram_legend.png", dpi=600)
