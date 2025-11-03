
# script to analyse bacterial growth rates
import polars as pl
import matplotlib.pyplot as plt
import seaborn.objects as so
from scipy.stats import mannwhitneyu, false_discovery_control, spearmanr
import warnings

# ~~~~~~~~~~~~ CONFIGS ~~~~~~~~~~~~ #
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "font.family": "Helvetica"
})

# ~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *


# ~~~~~~~~~~~~ MAPPING ~~~~~~~~~~~~ #
mapping = load_mapping()


# ~~~~~~~~~~~~ AESTHETICS ~~~~~~~~~~~~ #
color_df = load_colors()
colors_family = color_dict_family()
colors_strain = color_dict_strain()
markers = {"native": "s", "non-native": "D"}
linestyles = {"native": "-", "non-native": "--"}
title_dict = {
            'at_col0': "$At$-SC$^{At}$", 
            'lj_col0': "$Lj$-SC$^{At}$", 
            'at_gifu': "$At$-SC$^{Lj}$", 
            'lj_gifu': "$Lj$-SC$^{Lj}$", 
            }

# load growth rate table
growth_data = pl.read_csv("raw_data/insitu.tsv", separator="\t")


# first, plot the growth rates over time
growth_evolution = (growth_data

                        # filter for continous experiment
                        .filter(pl.col("sample_type") == "generation")

                        # add the color dataframe
                        .join(color_df, how="left", on="strain")

                        #sort, for the funsies
                        .sort("generation", "family", "plant", "syncom")
                        )


# normalize growth rates with respect to generation 1
norm_val = (growth_evolution
                .filter(pl.col("generation") == 1)
                .group_by("strain", "plant", "syncom").agg(pl.col("log2(PTR)").mean().alias("norm_val"))
                )


# perform statistics
plant_ls = []
syncom_ls = []
strain_ls = []
pval_ls = []
statistic_ls = []
for plant in ["col0", "gifu"]:
    for syncom in ["at", "lj"]:
        condi_subset = (growth_evolution
                            .filter(pl.col("plant") == plant)
                            .filter(pl.col("syncom") == syncom)
                            .filter(pl.col("generation") != 0)
                            .join(norm_val, how="left", on=["strain", "plant", "syncom"])
                            .with_columns((pl.col("log2(PTR)")/pl.col("norm_val")).alias("log2(PTR)"))
                            .filter(pl.col("log2(PTR)").is_not_null())
                        )
        for strain in condi_subset["strain"].unique():
            strain_susbet = condi_subset.filter(pl.col("strain") == strain)
            corr_res = spearmanr(strain_susbet["generation"].to_numpy(), 
                                 strain_susbet["log2(PTR)"].to_numpy())
            
            comp_res = mannwhitneyu(strain_susbet.filter(pl.col("generation") == pl.col("generation").min())["log2(PTR)"],
                                    strain_susbet.filter(pl.col("generation") == pl.col("generation").max())["log2(PTR)"])
            
            # gather data
            stats = comp_res
            plant_ls.append(plant)
            syncom_ls.append(syncom)
            strain_ls.append(strain)
            pval_ls.append(stats.pvalue)
            statistic_ls.append(stats.statistic)

stats_res = (pl.DataFrame({"plant": plant_ls, 
                           "syncom": syncom_ls, 
                           "strain": strain_ls, 
                           "pval": pval_ls,
                           "statistic": statistic_ls
                           })
                
                .filter(pl.col("pval").is_not_nan())
                .with_columns(pl.col("pval").map_batches(lambda x: false_discovery_control(x), return_dtype=pl.Float32).alias("padj"))
                .filter(pl.col("padj")<0.05)
            )


# once, plot the raw growth rates to highlight strain specific differences
print(f"n = {growth_evolution.filter(pl.col('generation')!=0).shape[0]}")
fig, axes = plt.subplots(nrows=4, figsize=(8,6))
idx=0
for plant in ["col0", "gifu"]:
    for syncom in ["at", "lj"]:
        condi_subset = (growth_evolution
                            .filter(pl.col("plant") == plant)
                            .filter(pl.col("syncom") == syncom)
                            .filter(pl.col("generation") != 0)
                            .join(norm_val, how="left", on=["strain", "plant", "syncom"])
                            .with_columns((pl.col("log2(PTR)")/pl.col("norm_val")).alias("log2(PTR)"))
                            .filter(pl.col("log2(PTR)").is_not_null())
                        )
        
        # get significant subset
        sig = stats_res.select("plant", "syncom", "strain").filter(pl.col("plant") == plant).filter(pl.col("syncom") == syncom).join(condi_subset, how="left", on = ["plant", "syncom", "strain"])
        

        ax = axes[idx]
        (so.Plot(data=condi_subset.to_pandas(),
                 x="generation",
                 y="log2(PTR)",
                 color="strain",
                 marker="treatment",)

                 # add artists
                 .add(so.Line(edgewidth=0, alpha=0.5), so.Agg(func="mean"), legend=False)
                 .add(so.Range(alpha=0.5), so.Est(func="mean", errorbar=("ci", 95)), legend=False)
                 .scale(color=colors_strain)
                 .scale(marker=markers)
                 .on(ax)
         
         ).plot()



        # general aesthetics
        if idx == 3:
            ax.set_xlabel("Plant cycle, $t$")
        else:
            ax.set_xlabel("")
        ax.set_ylabel("Growth rate, $r$")
        ax.set_xticks(range(1,17))
        ax.set_xlim((0.5, 16.5))
        ax.spines[["top", "right"]].set_visible(False)
        condition = title_dict[f'{syncom}_{plant}']
        ax.set_title(f"$In-situ$ growth rates across time for {condition}")
        ax.axhline(y=1, color="grey", alpha=0.5, linestyle="dashed")
        idx+=1

# save figure
fig.tight_layout(pad=0, h_pad=1)
fig.savefig(f"figures/shotgun_growth_over_time.pdf")
fig.savefig(f"figures/shotgun_growth_over_time.png", dpi=600)


