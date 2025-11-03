


# imports
import polars as pl
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, false_discovery_control
import seaborn.objects as so
from matplotlib.ticker import PercentFormatter
import warnings

# configs
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "font.family": "Helvetica"
})



from amplicon2_analysis_functions import *

# load the data
asv_table = pl.read_csv("raw_data/amplicon2_comp_otu.tsv", separator="\t")
design_table = pl.read_csv("raw_data/amplicon2_comp_design.tsv", separator="\t")
sc_strains = pl.read_csv("mapping/syncom_strains.tsv", separator="\t")
taxonomy = pl.read_csv("mapping/taxonomy.tsv", separator="\t")
evolved = pl.read_csv("mapping/evolved_strains.tsv", separator="\t")

# calculate RA
asv_table_RA = asv_table.with_columns(pl.col(sample) / pl.col(sample).sum() for sample in design_table["sampleID"])


# add meta data and melt asv table
asv_melted = (asv_table_RA
                        .unpivot(index="#Strain", on=design_table["sampleID"].to_list(), variable_name="sampleID", value_name="RA")
                        .join(design_table, how="left", on="sampleID"))




# define colors and stuff
color_dict = {"col0": "red", "gifu": "blue"}
order = design_table["treatment"].unique(maintain_order=True)
title_dict = {
            'at0_vs_lj0-at':       "$At$-SC$^{*}$\nvs\n$Lj$-SC$^{*}$", 
            'at(gifu)_vs_lj0-at':  "$At$-SC$^{Lj}$\nvs\n$Lj$-SC$^{*}$", 
            'at(col0)_vs_lj0-at':  "$At$-SC$^{At}$\nvs\n$Lj$-SC$^{*}$", 
            'at0_vs_lj(col0)-lj':  "$Lj$-SC$^{At}$\nvs\n$At$-SC$^{*}$", 
            'at0_vs_lj(gifu)-lj':  "$Lj$-SC$^{Lj}$\nvs\n$At$-SC$^{*}$", 
            'at0_vs_lj0-lj':    "$Lj$-SC$^{*}$\nvs\n$At$-SC$^{*}$", 
            }
order_treatment = [ 'at0_vs_lj0',
                    'at(gifu)_vs_lj0',
                    'at(col0)_vs_lj0',
                    'at0_vs_lj(col0)',
                    'at0_vs_lj(gifu)',]

xlab_dict = {"col0": "$At$ (Col-0)", "gifu": "$Lj$ (Gifu)"}
marker_dict = {'at0_vs_lj0': "o",
               'at(gifu)_vs_lj0':"D",
               'at(col0)_vs_lj0':"s",
               'at0_vs_lj(col0)':"D",
               'at0_vs_lj(gifu)': "s",}

# pad zeros 
min_val = asv_melted.filter(pl.col("treatment") != "input").filter(pl.col("RA")!= 0)["RA"].min()
max_val = asv_melted.filter(pl.col("treatment") != "input").filter(pl.col("RA")!= 0)["RA"].max()
padding_val = (min_val/max_val)*min_val



# now calculate host preference
host_pref_df = (asv_melted
                
                    # fill zeros
                    .with_columns(pl.col("RA").replace(0, padding_val))

                    # add syncom information for calculation of host preference
                    .with_columns(pl.when(pl.col("#Strain").str.starts_with("Root")).then(pl.lit("at")).otherwise(pl.lit("lj")).alias("strains"))
                    
                    # remove inputs
                    .filter(pl.col("treatment") != "input")

                    # calculate host preference
                    .group_by("sampleID", "plant", "strains", "treatment")
                    .agg(pl.col("RA").sum())
                    
                    # sort results
                    .sort("sampleID")
                    )



# now perform the plotting
for syncom in host_pref_df["strains"].unique():

    # determine native and non-native host
    native_host = "col0" if syncom == "at" else "gifu"
    non_native_host = "col0" if syncom == "lj" else "gifu"


    # perform statistics here
    pval_ls = []
    host_pref_ls = []
    for plant in [native_host, non_native_host]:

        # get mean to renormalize data
        ancestral = (host_pref_df
                                .filter(pl.col("strains")==syncom)
                                .filter(pl.col("plant") == plant)
                                .filter(pl.col("treatment") == "at0_vs_lj0")["RA"].median()
                                )
        
        stats_subest = (host_pref_df
                            .filter(pl.col("strains") == syncom)
                            .filter(pl.col("plant") == plant)
                            .with_columns(pl.col("RA")/ancestral)
                            )

        pvalue = mannwhitneyu(x=stats_subest.filter(pl.col("treatment") == "at0_vs_lj0")["RA"],
                              y=stats_subest.filter(pl.col("treatment").str.contains_any([f"{syncom}({plant})"]))["RA"],
                            ).pvalue
        pval_ls.append(pvalue)
        host_pref_ls.append(stats_subest)
    
    # adjust pvals and store in dict
    padj_dict = {key: val for key, val in zip((native_host, non_native_host), false_discovery_control(pval_ls))}
    plot_df = pl.concat(host_pref_ls)

    # create figure here
    fig, axes = plt.subplots(figsize=(4,3), ncols=2)

    for idx_plant, plant in enumerate(["col0", "gifu"]):

        # define correct axes
        ax=axes[idx_plant]

        # subset the data
        subset = (plot_df

                    # subsetting
                    .filter(pl.col("strains") == syncom)
                    .filter(pl.col("plant") == plant)
                    .filter(pl.col("treatment").str.contains_any(["at0_vs_lj0", f"{syncom}({plant})"]))
                    .sort(pl.col("treatment").cast(pl.Enum(order)))
                    )

        # fancy plotting looool
        print(f"syncom: {syncom}; plant: {plant} - n = {subset.shape[0]}")
        (so.Plot(data=subset.to_pandas(),
                    x="treatment",
                    y="RA",
                    color="plant",)

                    # add the individual data-points
                    .add(so.Dash(color="black", width=0.4), so.Agg(func="median"))
                    .add(so.Range(color="black"), so.Est(func="median", errorbar="sd"))
                    .add(so.Dot(alpha=0.6, pointsize=4, edgewidth=0), so.Jitter(width=0.25))

                    # adjust the elements
                    .scale(color=so.Nominal(color_dict))

                    # put on figure
                    .on(ax)
                    .plot()
        )
        ymax=subset["RA"].max()
        y_off=(subset["RA"].max()-ax.get_ylim()[0])*0.2
        pct_increase = round(((subset.filter(pl.col("treatment")!="at0_vs_lj0")["RA"].median() - subset.filter(pl.col("treatment")=="at0_vs_lj0")["RA"].median()) / subset.filter(pl.col("treatment")=="at0_vs_lj0")["RA"].median()) * 100, 2)

        # general figure aesthetics
        plant_str = "At" if plant == "col0" else "Lj"
        ax.set_ylabel("Bacterial fitness")
        ax.set_xlabel("")
        ax.set_title(f"Fitness on ${plant_str}$")
        ax.set_ylim((ax.get_ylim()[0], ymax+y_off))
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        ax.spines[["top", "right"]].set_visible(False)
        ax.axhline(y=1, linestyle="--", color="k", alpha=0.25)
        ax.set_xticklabels(["Ancestral", f"${syncom.title()}$-SC" + "$^{" + plant_str + "}$"])
        ax.text(x=0.5, 
                y=ymax, 
                s=f"{convert_pvalue_to_asterisks(padj_dict[plant])}\n+{pct_increase}%",
                ha="center")

        if idx_plant == 1:
            ax.set_ylabel("")

    # remove all legends
    [legend.set_visible(False) for legend in fig.legends]

    # title
    fig.suptitle(f"Fitness of evolved ${syncom.title()}$-SC")

    # save the figure
    fig.tight_layout()
    fig.savefig(f"figures/amplicon2_fitness_vs_ancestral_{syncom}.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(f"figures/amplicon2_fitness_vs_ancestral_{syncom}.png", bbox_inches="tight", pad_inches=0.05, dpi=600)

