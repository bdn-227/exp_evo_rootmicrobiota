
# script to replicate figure 2b from 
# Lenksi LTEE shotgun paper but on strain level
# https://www.nature.com/articles/nature24287/figures/2


# import modules
import polars as pl
import matplotlib.pyplot as plt
import numpy as np
import seaborn.objects as so
from math import ceil
from matplotlib import ticker
import warnings

# ~~~~~~~~~~~~ CONFIGS ~~~~~~~~~~~~ #
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "font.family": "Helvetica"
})



# ~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *


# ~~~~~~~~~~~~ aesthetics ~~~~~~~~~~~~ #
markers = {"native": "s", "non-native": "D"}
linestyles = {"native": "-", "non-native": ":"}


# ~~~~~~~~~~~~ MAPPING ~~~~~~~~~~~~ #
mapping = load_mapping()


# ~~~~~~~~~~~~ CREATE DUMMY WITH ALL ~~~~~~~~~~~~ #
full_mapping = pl.read_csv("mapping/shotgun_dummy.tsv", separator="\t").filter(pl.col("sample_type") == "generation").filter(pl.col("treatment")!="input")
dummy_df = (pl.concat([full_mapping.filter(pl.col("syncom") == syncom).join((pl.read_csv(f"mapping/syncom_mapping_{syncom}.txt").with_columns(pl.lit(syncom).alias("syncom")).rename({syncom: "strain"}).join(load_colors().select("strain", "family", "color").unique(), how="left", on="strain")), how="left", on="syncom") for syncom in ["at", "lj"]]))


# load the actual variant calls
variants_long = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv", 
                             separator="\t",
                             schema_overrides={"line": str}
                            )
                            .filter(pl.col("type")!="upstream")
                )



# genome sizes
genome_sizes = (pl.read_csv('raw_data/coverage_strainwise.tsv', separator="\t")
                . with_columns( ("contig_" + pl.col("#rname").str.split_exact("_", n=2).struct[-1]).alias("contig") )
                .select("strain", "contig", "endpos")
                .unique()
                .group_by("strain")
                .agg(pl.col("endpos").sum().alias("genomesize"))
                )



# same plot but on strain level
df_cols = ["line", "generation", "treatment", "strain", "ID", "AF"]
group_cols = ["line", "generation", "treatment", "strain"]
mt_df = (variants_long
            .select(df_cols)
            .group_by(group_cols)
            .agg(pl.col("AF").sum().alias("mt"))
            .sort("line", "strain", "generation")
            )

fixed = (variants_long
            .select(df_cols) 
            .filter(pl.col("AF") >= 0.95)
            .group_by(group_cols)
            .agg(pl.col("AF").sum().alias("fixed"))
            .sort("line", "generation")
            )


mt_df = (mt_df
            .join(fixed, how="left", on=group_cols)
            .fill_null(0)
            )

mt_df = dummy_df.join(mt_df, how="left", on=group_cols).with_columns((pl.col("family") + pl.col("line")).alias("identifier")).fill_null(0)



colors = {id: color for id, color in zip(mt_df["identifier"], mt_df["color"])}


# define the order
strains = (mt_df
                .with_columns(pl.when(pl.col("strain").str.starts_with("Root"))
                                .then(pl.lit("at"))
                                .otherwise(pl.lit("lj")).alias("syncom"))
                                .sort("family", "syncom")["strain"]
                                .unique(maintain_order=True)
                                )


# define figure layout
scaling_factor=0.8
figsize=(18*scaling_factor, 20*scaling_factor)
ncols=4
nrows = ceil(strains.len()/ncols)
fig, axes = plt.subplots(figsize=figsize, ncols=ncols, nrows=nrows)
dummy = np.arange(ncols*nrows).reshape(nrows, ncols)
xtick_idx = dummy.reshape(-1)[0:strains.len()][::-1][0:4]
ytick_idx = dummy[:, 0]
axes = axes.reshape(-1)


# get y limits
xlim = (0.25,16.25)
xticks = np.arange(1,17)
ylim = (0, int(mt_df["mt"].max()*1.1))
yticks = np.arange(0, np.max(ylim), 10)
print(f"mt for all strains: n = {mt_df.shape[0]}")
for strain_idx, strain in enumerate(strains):
    
    # define axis
    ax=axes[strain_idx]

    for treatment in mt_df["treatment"].unique():
        (so.Plot(data=mt_df.filter(pl.col("treatment") == treatment).filter(pl.col("strain")==strain).to_pandas(),
                x="generation",
                y="mt",
                color="identifier",
                marker="treatment",)
                
                # add the artists
                .add(so.Line(linestyle=linestyles[treatment], alpha=0.5, edgewidth=0))

                # scale artists
                .scale(color=colors)
                .scale(marker=markers)

                .on(ax)
                ).plot()


    # add labels
    if strain_idx in ytick_idx:
        ax.set_ylabel("Mutations, $M(t)$")
    else:
        ax.set_ylabel("")
    

    if strain_idx in xtick_idx:
        ax.set_xlabel("Plant cycle, $t$")
    else:
        ax.set_xlabel("")
    

    # reformat axis
    ax.set_xlim(xlim)
    ax.set_xticks(xticks)
    ax.set_ylim(ylim)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=10, integer=True, steps=[1,2,5,10]))

    # aesthetics
    ax.spines[["top", "right"]].set_visible(False)

    # title
    family = mt_df.filter(pl.col("strain")==strain)["family"][0]
    ax.set_title(f"{strain}\n(${family}$)")


# manage legend
[legend.set_visible(False) for legend in fig.legends]

# remove empty plots
[axes[e].set_visible(False) for e in np.arange(strain_idx+1, axes.shape[0])]


# write
fig.tight_layout(pad=1)
fig.savefig(f"figures/shotgun_mutation_mt_strain.pdf")
fig.savefig(f"figures/shotgun_mutation_mt_strain.png", dpi=600)
plt.clf()




# investigate Root142 outlier --> line 11
outlier_vals = (variants_long
                    .filter(pl.col("line")=="11")
                    .filter(pl.col("af_type")=="original")
                    .filter(pl.col("strain")=="Root142")
                    .filter(pl.col("generation").is_in([16])))
outlier_vals.group_by("desc").len()
outlier_vals.filter(pl.col("desc") == "adenylate cyclase [EC:4.6.1.1]")
# cardiolipin semms to regulate two component systems: https://journals.asm.org/doi/full/10.1128/iai.00046-23







# define figure layout
fig, axes = plt.subplots(figsize=figsize, ncols=ncols, nrows=nrows)
axes = axes.reshape(-1)


# get y limits
xlim = (-1, int(mt_df["mt"].max()*1.1))
xticks = np.arange(0, np.max(xlim), 10)
ylim = (-1, int(mt_df["fixed"].max()*1.1))

for strain_idx, strain in enumerate(strains):
    
    # define axis
    ax=axes[strain_idx]

    for treatment in mt_df["treatment"].unique():
        (so.Plot(data=mt_df.filter(pl.col("treatment") == treatment).filter(pl.col("strain")==strain).to_pandas(),
                 x="mt",
                 y="fixed",
                 color="identifier",
                 marker="treatment",)
                
                 # add the artists
                 .add(so.Line(linestyle=linestyles[treatment], 
                              alpha=0.5, 
                              edgewidth=0, 
                              pointsize=3))

                  # scale artists
                  .scale(color=colors)
                  .scale(marker=markers)

                  .on(ax)
                  ).plot()
        
        # add labels
        if strain_idx in ytick_idx:
            ax.set_ylabel("Fixed mutations")
        else:
            ax.set_ylabel("")
        

        if strain_idx in xtick_idx:
            ax.set_xlabel("Mutations, $M(t)$")
        else:
            ax.set_xlabel("")

        # aesthetics
        ax.spines[["top", "right"]].set_visible(False)

        # title
        family = mt_df.filter(pl.col("strain")==strain)["family"][0]
        ax.set_title(f"{strain}\n(${family}$)")



# manage legend
[legend.set_visible(False) for legend in fig.legends]

# remove empty plots
[axes[e].set_visible(False) for e in np.arange(strain_idx+1, axes.shape[0])]


# write
fig.tight_layout(pad=1)
fig.savefig(f"figures/shotgun_mutation_mt_vs_fixed_strain.pdf")
fig.savefig(f"figures/shotgun_mutation_mt_vs_fixed_strain.png", dpi=600)
plt.clf()

