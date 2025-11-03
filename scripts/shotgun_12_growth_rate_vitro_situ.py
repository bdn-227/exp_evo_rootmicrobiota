
# script to analyse bacterial growth rates
# specically, we are investigating, whether bacterial growth rates
# (in situ and in vitro) are correlated with ne substitution rate


import polars as pl
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import seaborn.objects as so
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import BaggingRegressor
import warnings
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "font.family": "Helvetica"
})

# ~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *


# ~~~~~~~~~~~~ MAPPING ~~~~~~~~~~~~ #
mapping = load_mapping()
color_df = load_colors()
colors_family = color_dict_family()
colors_strain = color_dict_strain()



# do the same with the invitro growth data
# (that is growth cruve from cultivation)
invitro_growth = pl.read_csv("raw_data/invitro.tsv", separator="\t")


# test for correlation between in-situ and in-vitro growth
merged_growth = (invitro_growth
                        .join(pl.read_csv("raw_data/insitu.tsv", separator="\t").group_by("strain").agg(pl.col("log2(PTR)").mean()), how="left", on="strain")
                        .filter(pl.col("log2(PTR)").is_not_null())    
                )

corr_res = spearmanr(merged_growth["r"].to_numpy(), 
                     merged_growth["log2(PTR)"].to_numpy())
corr_res


X = merged_growth['r'].to_numpy().reshape(-1,1)
y = merged_growth["log2(PTR)"].to_numpy()
model = BaggingRegressor(LinearRegression(),
                         n_estimators=100,
                         max_samples=1.0,
                         bootstrap=True)

model.fit(X, y)
bootstrapped_preds = (pl.DataFrame([m.predict(X) for m in model.estimators_])
                            .with_columns(pl.lit(X.reshape(-1))
                                          .alias("r"))
                                          .unpivot(index="r")
                                          )


# plot scatter
print(f"n = {merged_growth.shape[0]}")
fig, axes = plt.subplots(figsize=(4,3))
(so.Plot(data=merged_growth.to_pandas(),
         x="r",
         y="log2(PTR)",
         )
         .add(so.Dot(alpha=0.5, edgewidth=0), color="strain", legend=False)
         .scale(color=colors_strain)
         .on(axes)
    
).plot()

(so.Plot(data=bootstrapped_preds.to_pandas(),
         x="r",
         y="value",
        )

        .add(so.Band(), so.Est(errorbar=('sd',2), n_boot=1000))
        .add(so.Line(color="k"), so.PolyFit(order=1))
        
        .on(axes)
        ).plot()


# general aesthetics
axes.set_yscale('log')
axes.set_xscale('log')
axes.spines[["top", "right"]].set_visible(False)
axes.set_title(f"Correlation analysis of growth rates")
axes.set_xlabel("Growth rate, $r$ (in-vitro)")
axes.set_ylabel("Growth rate, $r$ (in-situ)")

# define axis boundries
xlim = axes.get_xlim()
ylim = axes.get_ylim()
ylim_min = np.min([xlim, ylim])
ylim_max = np.max([xlim, ylim])

# write the statistics in condition specific manner
axes.text(x=merged_growth["r"].min(),
          y=np.max(axes.get_ylim()) * 0.8,
          s=f"$r²={round(corr_res.statistic, 3)}$; $p=${pretty_exponent(corr_res.pvalue, 2)}\n$n$ = {merged_growth.shape[0]}",
          ha="left",
          va="top")


    
fig.tight_layout(pad=0)
fig.savefig(f"figures/shotgun_growth_rate_correlation_vitro_situ.pdf")
fig.savefig(f"figures/shotgun_growth_rate_correlation_vitro_situ.png", dpi=600)


