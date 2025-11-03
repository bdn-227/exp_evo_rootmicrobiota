
# script to plot the host preference as fitness gain from
# generation 1 to 16 

import polars as pl
import seaborn.objects as so
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import false_discovery_control, mannwhitneyu
from matplotlib.patches import Patch
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


# styles
color_dict = {"col0": "red", "gifu": "blue"}
marker_dict = {"ancestral": "o", "native": "s", "non_native": "D"}
title_dict = {"col0": "$A. thaliana$ (Col-0)", "gifu": "$L. japonicus$ (Gifu)"}
legend_color = [Patch(label="Fitness on", alpha=0)] + [Patch(facecolor=color, label=title_dict[label], alpha=1) for label, color in color_dict.items()]


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
                    .with_columns(pl.col("RA")/pl.col("depth"))
                    .group_by("sampleID", "plant", "strains", "gen", "condition")
                    .agg(pl.col("RA").sum())
                    .with_columns(pl.lit("miseq").alias("dataset"))
                    .sort("sampleID")
                    )


# normalize the data with respect to native host
syncom_ls = []
timepoint_ls = []
condition_ls = []
pval_ls = []
df_list = []
ymax_ls=[]
for syncom in ["at", "lj"]:
    for timepoint in host_pref_df["gen"].unique():

        # subset to get the conditions
        host_pref_subset = host_pref_df.filter(pl.col("strains") == syncom).filter(pl.col("gen") == timepoint)
        
        for condition in ["native", "non_native"]:
            
            # subset the data
            other = "native" if condition == "non_native" else "non_native"
            condition_subset = host_pref_subset.filter(pl.col("condition") != other)


            # now, re-normalize the data with respect to evolving host
            if condition == "non_native":
                non_evolving_host = "col0" if syncom == "at" else "gifu"
            else:
                non_evolving_host = "col0" if syncom == "lj" else "gifu"
            
            condition_subset = (condition_subset
                                        .with_columns(pl.lit(condition_subset.filter(pl.col("plant") == non_evolving_host)["RA"].mean()).alias("meanRA"))
                                        .with_columns( (pl.col("RA") / pl.col("meanRA")).alias("normalizedRA") )
                                        )
            
            # now perform statistics
            native_set = condition_subset.filter(pl.col("plant") == non_evolving_host)["normalizedRA"].to_numpy()
            nnative_set= condition_subset.filter(pl.col("plant") != non_evolving_host)["normalizedRA"].to_numpy()
            pval = mannwhitneyu(native_set, nnative_set).pvalue
            maxval = np.array([native_set.max(), nnative_set.max(), np.std(native_set) + np.mean(native_set), np.std(nnative_set) + np.median(nnative_set)])
            ymax = maxval.max()+0.1
            
            # store results
            syncom_ls.append(syncom)
            timepoint_ls.append(timepoint)
            condition_ls.append(condition)
            pval_ls.append(pval)
            df_list.append(condition_subset)
            ymax_ls.append(ymax)


# adjust pvalues
padj = false_discovery_control(pval_ls)
sig = [convert_pvalue_to_asterisks(pval) for pval in pval_ls]

# concat data together
host_pref_df = pl.concat(df_list)
stats_df = (pl.DataFrame({"strains": syncom_ls,
                          "gen": timepoint_ls,
                          "condition": condition_ls,
                          "pval": pval_ls,
                          "ymax": ymax_ls})
                          .with_columns(pl.Series(name="padj", values=padj))
                          .with_columns(pl.Series(name="significant", values=sig))
                          )



# now perform the plotting
for syncom in ["at", "lj"]:

    # get host
    nnative_host = "col0" if syncom == "lj" else "gifu"

    # initialize the figure
    sc_subset = host_pref_df.filter(pl.col("strains") == syncom)
    fig, axes = plt.subplots(
                             ncols=2,
                             nrows=sc_subset["dataset"].unique().len(), 
                             figsize=(10,3))
    
    # get ylims
    ylim=(sc_subset["normalizedRA"].min()*.8, sc_subset["normalizedRA"].max()*1.2)


    for idx_dataset, dataset in enumerate(sc_subset["dataset"].unique()):

        # subset the data
        dataset_subset = (sc_subset
                            .filter(pl.col("dataset") == dataset)
                            .sort("condition", descending=True)
                            )


        for idx_condition, condition in enumerate(["non_native", "native"]):
            
            # subset the data
            other = "native" if condition == "non_native" else "non_native"
            condition_subset = (dataset_subset
                                    .filter(pl.col("condition") != other)
                                    .with_columns(pl.lit(condition).alias("condition"))
                                    .sort("gen")
                                    )

            # define the correct axis
            ax = axes.reshape((-1 ,1))[idx_condition, idx_dataset]

            # plot
            print(f"syncom: {syncom} - {dataset} - {condition} - n = {condition_subset.shape[0]}")
            (so.Plot(data=condition_subset.to_pandas(),
                        x="gen",
                        y="normalizedRA",
                        color="plant",
                        marker="condition")

                        # add the plot elements
                        .add(so.Dot(alpha=0.3, pointsize=4, marker=marker_dict[condition], edgewidth=0), so.Jitter(), so.Dodge(), legend=False)
                        .add(so.Range(), so.Est(errorbar="sd"), so.Dodge(), legend=False)
                        .add(so.Dot(marker=marker_dict[condition]), so.Agg("median"), so.Dodge(), legend=False)


                        # adjust the elements
                        .scale(color=so.Nominal(color_dict))

                        # put on figure
                        .on(ax)
                        .plot()
            )

            # calculate fitness gain here
            if condition == "non_native":
                evolving_host = "col0" if syncom == "lj" else "gifu"
                non_evolving_host = "col0" if syncom == "at" else "gifu"
            else:
                evolving_host = "col0" if syncom == "at" else "gifu"
                non_evolving_host = "col0" if syncom == "lj" else "gifu"

            # add the fitted line
            initial = condition_subset.filter(pl.col("plant") == evolving_host).filter(pl.col("gen")==1)["normalizedRA"].median()
            final   = condition_subset.filter(pl.col("plant") == evolving_host).filter(pl.col("gen")==16)["normalizedRA"].median()
            pct_change = round(((final-initial)/initial)*100, 1)

            # add the line connecting
            ax.axline((1, initial), 
                      (16, final), 
                      #linewidth=4, 
                      color=color_dict[evolving_host])
            ax.spines[["top", "right"]].set_visible(False)


            # add statistics
            stats_subset = (stats_df
                            .filter(pl.col("strains") == syncom)
                            .filter(pl.col("condition") == condition)
                            )
            
            # make title
            ax.set_ylim((ax.get_ylim()[0], stats_subset["ymax"].max()+0.1))
            ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
            ax.set_title(f"Fitness of {condition}ly evolved ${syncom.title()}$-SC (${syncom.title()}$-SC$^" + "{" + ("At" if evolving_host == "col0" else "Lj") + "}$)")
            ax.set_ylabel("Bacterial fitness")
            ax.set_xticks(condition_subset["gen"].unique())
            ax.axhline(y=1, 
                       color="k", 
                       linestyle="-.",
                       alpha=0.5)
            
            # manage xlabel
            if condition == "native":
                ax.set_xlabel("Plant cycle, $t$")
            else:
                ax.set_xlabel("")
            
            for row in stats_subset.iter_rows(named=True):
                ax.text(x=row["gen"],
                        y=row["ymax"],
                        s=row["significant"],
                        ha="center",
                        va="center")

    fig.tight_layout()
    fig.savefig(f"figures/amplicon_fitness_gain_all_{syncom}.pdf")
    fig.savefig(f"figures/amplicon_fitness_gain_all_{syncom}.png", dpi=600)




