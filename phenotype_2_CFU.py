

# script to compare the ratio of plant:bacterial reads
# to assess whether new ecological niches arose during the experiment



# ~~~~~~~~~~~~ IMPORTS ~~~~~~~~~~~~ #
import polars as pl
import matplotlib.pyplot as plt
import seaborn.objects as so
from scipy.stats import spearmanr, wilcoxon, false_discovery_control
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import BaggingRegressor
import warnings

# configs
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "font.family": "Helvetica"
})


# ~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *

# aesthetics
markers={"native": "s", "non-native": "D"}
colors = {
    "1": "darkred",
    "2": "red",
    "3": "orangered",
    "4": "coral",
    "5": "salmon",
    "6": "black",
    "7": "dimgrey",
    "8": "slategrey",
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





# same type of analysis for CFUs obtained from spectral scans
reisolations = (pl.read_csv("raw_data/reisolation.tsv", separator="\t", schema_overrides={"line": str}))

# write table with the dilution rates
dilution_rates = (reisolations
                    .select("filename", "generation", "line", "replicate", "date", "dilution")
                    .unique()
                    .with_columns(pl.col("line").cast(float).alias("line_f"))
                    .sort("generation", "line_f", "replicate")
                    .drop("line_f")
                    )
dilution_rates.write_csv("results/S_dilution_rates.tsv", separator="\t")


# now, determine the CFUs per ml root isolate
cfu = (reisolations
                .filter(~pl.col("strainID").is_in(["medium", "blank"]))
                .group_by("line", "generation", "replicate", "syncom", "dilution", "treatment")
                .len("numstrains")
                .with_columns((pl.col("numstrains")*pl.col("dilution")).alias("CFU"))
                )



X = cfu['generation'].to_numpy().reshape(-1,1)
y = cfu['CFU'].to_numpy()
model = BaggingRegressor(LinearRegression(),
                         n_estimators=100,
                         max_samples=1.0, # 100% of the dataset
                         bootstrap=True)

model.fit(X, y)
bootstrapped_preds = (pl.DataFrame([m.predict(X) for m in model.estimators_])
                            .with_columns(pl.lit(X.reshape(-1))
                                          .alias("generation"))
                                          .unpivot(index="generation")
                                          )

# also test for correlation
corr_res = spearmanr(cfu["generation"].to_numpy(), 
                     cfu["CFU"].to_numpy())



# always test for significance against ancestral
pval_ls = []
gen_ls = []
gen_1 = cfu.filter(pl.col("generation") == 1).group_by("line").agg(pl.col("CFU").mean().log(base=10))
for gen in cfu["generation"].unique().sort()[1:]:
    gen_x = cfu.filter(pl.col("generation") == gen).group_by("line").agg(pl.col("CFU").mean().log(base=10))
    gen_x = gen_x.join(gen_1.rename({"CFU": "CFU_1"}), how="inner", on="line")
    comp_res = wilcoxon(gen_x["CFU_1"], gen_x["CFU"])
    pval_ls.append(comp_res.pvalue)
    gen_ls.append(gen)
stat_res = (pl.DataFrame({"generation": gen_ls, 
                          "pval": pval_ls,})
                    .with_columns(pl.col("pval").map_batches(lambda x: false_discovery_control(x)).alias("padj"))
                    .with_columns(pl.col("padj").map_elements(lambda x: convert_pvalue_to_asterisks(x), return_dtype=str).alias("sig"))
           )


# plot the absolute abundances
fig, axes = plt.subplots(figsize=(8,3))

# perform the plotting
print(f"cfu: n = {cfu.shape[0]}")
(so.Plot(data=cfu.to_pandas(),
         x="generation",
         y="CFU",
        )
        .add(so.Dot(edgewidth=0, alpha=0.5), so.Jitter(), color="line", marker="treatment", legend=False)
        .add(so.Line(color="k"), so.PolyFit(order=1))
        .scale(color=colors)
        .scale(marker=markers)
        .on(axes)
        ).plot()


(so.Plot(data=bootstrapped_preds.to_pandas(),
         x="generation",
         y="value",
        )
        .add(so.Band(), so.Est(errorbar=('sd',2), n_boot=1000))
        .on(axes)
        ).plot()


# manage axis
axes.set_xticks(range(1,17))
axes.set_yscale('log')
axes.set_ylim((axes.get_ylim()[0], axes.get_ylim()[1]*1.5))

# labels
axes.set_xlabel("Plant cycle, $t$")
axes.set_ylabel("CFU")
axes.set_title("Bacterial load")

# add statistic
axes.text(x=1,
          y=axes.get_ylim()[1]*0.7,
          s=f"$r²={round(corr_res.statistic, 3)}$; $p=${pretty_exponent(corr_res.pvalue, 2)}",
          ha="left")


# general aesthetics
axes.spines[["top", "right"]].set_visible(False)

fig.tight_layout(pad=0.2)
fig.savefig(f"figures/abs_abundance_cfu.pdf")
fig.savefig(f"figures/abs_abundance_cfu.png")


