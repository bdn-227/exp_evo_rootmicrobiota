
# script to analyse bacterial growth rates
# specically, we are investigating, whether bacterial growth rates
# (in situ and in vitro) are correlated with substitution rate


import polars as pl
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import seaborn.objects as so
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import BaggingRegressor
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


# ~~~~~~~~~~~~ plotting aesthetics ~~~~~~~~~~~~ #
color_df = load_colors()
colors_family = color_dict_family()
colors_strain = color_dict_strain()
markers = {"native": "s", "non-native": "D"}
colors_plant = {"col0": "red", "gifu": "blue"}
linestyles = {"native": "-", "non-native": "--"}
title_dict = {
            'at_col0': "$At$-SC$^{Col-0}$", 
            'lj_col0': "$Lj$-SC$^{Col-0}$", 
            'at_gifu': "$At$-SC$^{Gifu}$", 
            'lj_gifu': "$Lj$-SC$^{Gifu}$", 
            }

# load growth rate table
growth_data = pl.read_csv("raw_data/insitu.tsv", separator="\t")


# mean get mean growth rate over time, by condition
aggregate_across = ["plant", "syncom", "treatment", "strain", "line"]
growth_evolution = (growth_data

                        # filter for continous experiment
                        .filter(pl.col("sample_type") == "generation")
                        .filter(pl.col("generation") != 0)

                        # add the color dataframe
                        .join(color_df, how="left", on="strain")

                        # aggregate
                        .group_by(aggregate_across + ["family"])
                        .agg(pl.col("log2(PTR)").mean())

                        )


# and import the inputs
variants_long = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv", 
                            separator="\t",
                            schema_overrides={"line": str}
                            )
                            .filter(pl.col("type")!="upstream")
                            .filter(pl.col("AF") > 0)
                            .with_columns(("g" + pl.col("generation").cast(str) + "_l" + pl.col("line")).alias("sampleID"))
                )


# get number of mutations
variants_n = variants_long.group_by(aggregate_across + ["generation"]).len("numvar").filter(pl.col("generation")==pl.col("generation").max())
merged_data = (growth_evolution
                    .join(variants_n, how="left", on=aggregate_across, validate="1:1")
                    .filter(pl.col("numvar").is_not_null())
                    .with_columns((pl.col("syncom") + "_" + pl.col("plant")).alias("condi"))
                    )


# test for correlation
corr_res = spearmanr(merged_data["log2(PTR)"], 
                     merged_data["numvar"])



X = merged_data["log2(PTR)"].to_numpy().reshape(-1,1)
y = merged_data["numvar"].to_numpy()
model = BaggingRegressor(LinearRegression(),
                         n_estimators=100,
                         max_samples=1.0, # 100% of the dataset
                         bootstrap=True)

model.fit(X, y)
bootstrapped_preds = (pl.DataFrame([m.predict(X) for m in model.estimators_])
                            .with_columns(pl.lit(X.reshape(-1))
                                          .alias("log2(PTR)"))
                                          .unpivot(index="log2(PTR)")
                                          )


# make a plot
print(f"n = {merged_data.shape[0]}")
fig, axes = plt.subplots(figsize=(4,3))
(so.Plot(data=merged_data.to_pandas(),
         x="log2(PTR)",
         y="numvar",
         color="strain",
         )
         
         .add(so.Dot(),  marker="treatment", legend=False)
         .scale(color=colors_strain)
         .on(axes)
         ).plot()

(so.Plot(data=bootstrapped_preds.to_pandas(),
         x="log2(PTR)",
         y="value",
        )

        .add(so.Band(color="k"), so.Est(errorbar=('sd',2), n_boot=1000))
        .add(so.Line(color="k"), so.PolyFit(order=1), legend=False)
        
        .on(axes)
        ).plot()


# general aesthetics
axes.spines[["top", "right"]].set_visible(False)
axes.set_title(f"Correlation between $in$-$situ$ $r$ and $M(t)$")
axes.set_xlabel("Growth rate, $r$ ($in$-$situ$)")
axes.set_ylabel("Mutations, $m(t)$")
axes.set_ylim((-1, np.max(axes.get_ylim())*1.1))
axes.set_xscale('log')

# write the statistics in condition specific manner
axes.text(x=merged_data["log2(PTR)"].min(),
          y=np.max(axes.get_ylim()) * 0.95,
          s=f"$r²={round(corr_res.statistic, 3)}$; $p=${pretty_exponent(corr_res.pvalue, 2)}\n$n$ = {merged_data.shape[0]}",
          ha="left",
          va="top")



fig.tight_layout(pad=0)
fig.savefig(f"figures/shotgun_correlation_insitu_mutations.pdf")
fig.savefig(f"figures/shotgun_correlation_insitu_mutations.png", dpi=600)


