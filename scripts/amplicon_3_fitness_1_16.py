
# script to plot the host preference as fitness gain from
# generation 1 and 16 only

import polars as pl
import seaborn.objects as so
import matplotlib.pyplot as plt
from scipy.stats import false_discovery_control, mannwhitneyu
from matplotlib import ticker
from matplotlib.ticker import PercentFormatter


# define fonts
plt.rcParams.update({
    "font.family": "Helvetica"
})

# ~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~ #
def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "n.s."


def percent_increase(old, new):
    if old == 0:
        raise ValueError("Old value cannot be zero for percent increase calculation.")
    return ((new - old) / old) * 100

# ~~~~~~~~~~~~ LOAD DATA ~~~~~~~~~~~~ #
design = pl.read_csv("raw_data/amplicon_evo_design.tsv", separator="\t")
otu = pl.read_csv("raw_data/amplicon_evo_otu.tsv", separator="\t")

# calculate host-preference
otu_long = otu.unpivot(index="strain", value_name="RA", variable_name="sampleID")
host_pref_df = (otu_long
                    .join(otu_long.group_by("sampleID").agg(pl.col("RA").sum().alias("depth")), how="left", on="sampleID")
                    .with_columns( pl.when(pl.col("strain").str.starts_with("Root")).then(pl.lit("at")).when(pl.col("strain").str.starts_with("Lj")).then(pl.lit("lj")).alias("strains") )
                    .join(design.select("plant", "condition", "gen", "sampleID", "compartment", "sc"), how="left", on="sampleID")
                    .filter(pl.col("compartment") == "root")
                    .filter(pl.col("sc") == "mixed")
                    .filter(pl.col("gen").is_in([1, 16]))
                    .with_columns(pl.col("RA")/pl.col("depth"))
                    .group_by("sampleID", "plant", "strains", "gen", "condition")
                    .agg(pl.col("RA").sum())
                    .sort("sampleID")
                    )


# normalize the data with respect to native host
syncom_ls = []
timepoint_ls = []
condition_ls = []
pval_ls = []
df_list = []
for syncom in ["at", "lj"]:
    for timepoint in host_pref_df["gen"].unique():

        # subset to get the conditions
        host_pref_subset = host_pref_df.filter(pl.col("strains") == syncom).filter(pl.col("gen") == timepoint)
        
        for condition in host_pref_subset["condition"].unique():
            
            # subset the data
            condition_subset = host_pref_subset.filter(pl.col("condition") == condition)


            # now, re-normalize the data with respect to native host
            native_host = "col0" if syncom == "at" else "gifu"
            condition_subset = (condition_subset
                                        .with_columns(pl.lit(condition_subset.filter(pl.col("plant") == native_host)["RA"].median()).alias("medianRA"))
                                        .with_columns( (pl.col("RA") / pl.col("medianRA")).alias("normalizedRA") )
                                        )
            
            # now perform statistics
            native_set = condition_subset.filter(pl.col("plant") == native_host)["normalizedRA"].to_numpy()
            nnative_set= condition_subset.filter(pl.col("plant") != native_host)["normalizedRA"].to_numpy()
            pval = mannwhitneyu(native_set, nnative_set).pvalue
            
            # store results
            syncom_ls.append(syncom)
            timepoint_ls.append(timepoint)
            condition_ls.append(condition)
            pval_ls.append(pval)
            df_list.append(condition_subset)


# adjust pvalues
padj = false_discovery_control(pval_ls)
sig = [convert_pvalue_to_asterisks(pval) for pval in padj]

# concat data together
host_pref_df = pl.concat(df_list)
stats_df = (pl.DataFrame({"strains": syncom_ls,
                          "gen": timepoint_ls,
                          "condition": condition_ls,
                          "pval": pval_ls})
                          .with_columns(pl.Series(name="padj", values=padj))
                          .with_columns(pl.Series(name="significant", values=sig))
                          )



# styles
color_plant = {"col0": "red", "gifu": "blue"}
color_syncom = {"at": "red", "lj": "blue"}
marker_dict = {"native": "s", "non_native": "D"}

title_dict = {"non_native": "$At$-SC$^{Lj}$ vs $Lj$-SC$^{At}$\n(Non-Native)",
              "native": "$At$-SC$^{At}$ vs $Lj$-SC$^{Lj}$\n(Native)",
              }

plant_dict = {"col0": "$A. thaliana$ (Col-0)",
              "gifu": "$L. japonicus$ (Gifu)"}



# get the actual min and max
at_sc = (pl.col("strains") == "at") & (pl.col("plant") == "gifu")
lj_sc = (pl.col("strains") == "lj") & (pl.col("plant") == "col0")
ymin = host_pref_df.filter(at_sc|lj_sc)["normalizedRA"].min()*0.9
ymax = host_pref_df.filter(at_sc|lj_sc)["normalizedRA"].max()*1.1


for condition in ["non_native", "native"]:

    # subset the data
    treatment_subset = (host_pref_df
                            .filter(pl.col("condition") == condition)
                            .filter(at_sc|lj_sc)
                            .with_columns( (pl.col("strains") + "_" + pl.col("gen").cast(str)).alias("x") )
                            .sort(pl.col("strains"), pl.col("gen"), pl.col("plant").cast(pl.Enum(["col0", "gifu"])))
                       )

    # print n's
    print(f"n's in the competition experiments {condition}: ")
    print(treatment_subset.group_by("strains").len())

    # make a nested grid spec
    fig, axes = plt.subplots(figsize=(4,4), ncols=1, nrows=1)

    # now perform plotting
    (so.Plot(data=treatment_subset.to_pandas(),
             x="x",
             y="normalizedRA",
             color="plant",
             marker="condition")

             # add the individual data-points
             .add(so.Bar(alpha=0.5, edgewidth=0), so.Agg(func="median"), legend=False)
             .add(so.Range(color="black"), so.Est(func="median", errorbar="sd"), legend=False)
             .add(so.Dot(alpha=0.8, pointsize=4, marker=marker_dict[condition]), so.Jitter(width=0.2), legend=False)

             # adjust the elements
             .scale(color=so.Nominal(color_plant))

             # put on figure
             .on(axes)
             .plot()
    )

        
    # aesthetics
    axes.set_xlabel("")
    axes.set_xticklabels(["#1", "#16", "#1", "#16", ])
    
    axes.set_ylim((ymin, ymax))
    axes.yaxis.set_major_locator(ticker.MultipleLocator(.2))
    axes.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    axes.set_ylabel("Bacterial fitness on non-host")
    axes.set_title(title_dict[condition])
    axes.spines[["top", "right"]].set_visible(False)
    axes.axhline(y=1, 
                 color="k", 
                 linestyle="-.",
                 zorder=100,
                 alpha=0.6
                )
    sec = axes.secondary_xaxis(location=0)
    axes.get_xticks()
    sec.set_xticks([0.5, 2.5], labels=['\n$At$-SC', '\n$Lj$-SC'])
    sec.tick_params(axis='x', which="both", length=5, width=0)
    
    # lines between the classes:
    sec2 = axes.secondary_xaxis(location=0)
    sec2.set_xticks([-0.5, 1.5, 3.5], labels=[])
    sec2.tick_params(axis='x', which="major", length=30, width=1)

    # now add stats
    pval_ls = []
    for syncom in ["at", "lj"]:
        sc_subset = treatment_subset.filter(pl.col("strains") == syncom)
        x = sc_subset.filter(pl.col("gen") == 1)["normalizedRA"]
        y = sc_subset.filter(pl.col("gen") == 16)["normalizedRA"]
        pval_ls.append(mannwhitneyu(x,y).pvalue)

        # print the percent gain for text
        print(f"fitness increase {syncom}-{condition}: {round(percent_increase(x.median(), y.median()), 2)}%")

    sig = [convert_pvalue_to_asterisks(padj) for padj in false_discovery_control(pval_ls)]

    
    y_vals = treatment_subset.group_by("strains").agg(pl.col("normalizedRA").max()+0.1)
    for idx, syncom in enumerate(["at", "lj"]):
        s=sig[idx]
        idx=idx*2
        y_val = y_vals.filter(pl.col("strains") == syncom)["normalizedRA"].max()
        axes.plot((idx, idx+1), [y_val]*2, color="k")
        axes.plot((idx, idx), [y_val, y_val-0.025], color="k")
        axes.plot((idx+1, idx+1), [y_val, y_val-0.025], color="k")
        axes.text(x=idx+0.5, y=y_val+0.025, s=s, ha="center")


    # title
    fig.tight_layout(h_pad=0)
    fig.savefig(f"figures/amplicon_fitness_gain_1_16_{condition}.pdf", bbox_inches="tight", pad_inches=0)
    fig.savefig(f"figures/amplicon_fitness_gain_1_16_{condition}.png", dpi=600, bbox_inches="tight", pad_inches=0)





