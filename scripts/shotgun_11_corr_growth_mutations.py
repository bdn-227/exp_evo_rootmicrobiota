
# script to analyse bacterial growth rates
# specically, we are investigating, whether bacterial growth rates
# (in situ and in vitro) are correlated with ne substitution rate


import polars as pl
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import seaborn.objects as so
import warnings


# ~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "font.family": "Helvetica"
})

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


# and import the inputs
variants_long = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv", 
                             separator="\t",
                             schema_overrides={"line": str}
                            )
                            .filter(pl.col("type")!="upstream")
                            .filter(pl.col("AF") > 0)
                            .with_columns(("g" + pl.col("generation").cast(str) + "_l" + pl.col("line")).alias("sampleID"))
                )



# load abundance data and calculate fold change
design = pl.read_csv("raw_data/amplicon_evo_design.tsv", separator="\t")
otu = pl.read_csv("raw_data/amplicon_evo_otu.tsv", separator="\t")
abundance_data = (otu
                    .unpivot(index="strain", value_name="RA", variable_name="sampleID")
                    .join(design, how="left", on="sampleID")
                    .filter(pl.col("compartment")=="root")
                    .filter(pl.col("sc").is_in(["at", "lj"]))
                    )


# calcluate fold-change
line_ls = []
strain_ls = []
fc_ls = []
mt_ls = []
for line in abundance_data["line"].unique():
    for strain in abundance_data["strain"].unique():

        try:

            # extract variants
            variant_subset = (variants_long
                                    .filter(pl.col("line")==line)
                                    .filter(pl.col("strain")==strain)
                                    .filter(pl.col("generation") == pl.col("generation").max())
                            )
            mt = variant_subset.shape[0]
            
            # extract abundance data
            abundance_subset = (abundance_data
                                    .filter(pl.col("line")==line)
                                    .filter(pl.col("strain")==strain)
                                    .sort("gen")
                                    )
            fc = abundance_subset.filter(pl.col("gen")==pl.col("gen").max())["RA"].max()/abundance_subset.filter(pl.col("gen")==pl.col("gen").min())["RA"].max()
            
            # aggregate data
            line_ls.append(line)
            strain_ls.append(strain)
            mt_ls.append(mt)
            fc_ls.append(fc)
        except:
            pass


merged_data = (pl.DataFrame({"line": line_ls, 
                             "strain": strain_ls, 
                             "mt": mt_ls, 
                             "fc": fc_ls,
                             })
                    .filter(pl.col("mt")!=0)
                    .filter(pl.col("fc")!=0)
                    .with_columns(pl.col("fc").log(base=2))
                    .sort("mt")
              )


# test for correlation
corr_res = spearmanr(merged_data["fc"], 
                     merged_data["mt"])


# make a plot
print(f"n = {merged_data.shape[0]}")
fig, axes = plt.subplots(figsize=(4,3))
(so.Plot(data=merged_data.to_pandas(),
         x="fc",
         y="mt",
         color="strain",
         marker="line",
         )
         
         .add(so.Dot(), legend=False)
         .scale(color=colors_strain)
         .on(axes)
).plot()



# general aesthetics
axes.spines[["top", "right"]].set_visible(False)
axes.set_title(f"Correlation between Log₁₀FC(RA) and $M(t)$")
axes.set_xlabel("log₁₀(RA fold change)")
axes.set_ylabel("Mutations, $M(t)$")


# write the statistics in condition specific manner
axes.text(x=merged_data["fc"].min(),
          y=np.max(axes.get_ylim()) * 0.9,
          s=f"$r²={round(corr_res.statistic, 3)}$; $p=${pretty_exponent(corr_res.pvalue, 2)}\n$n$ = {merged_data.shape[0]}",
          ha="left",
          va="top")

# Plot regression line
xseq = np.linspace(merged_data["fc"].min(), merged_data["fc"].max(), num=100)
b, a = np.polyfit(merged_data["fc"], merged_data["mt"], deg=1)
axes.plot(xseq, a + b * xseq, color="k", lw=1)


fig.tight_layout(pad=0)
fig.savefig(f"figures/shotgun_correlation_abundance_mutations.pdf")
fig.savefig(f"figures/shotgun_correlation_abundance_mutations.png", dpi=600)


