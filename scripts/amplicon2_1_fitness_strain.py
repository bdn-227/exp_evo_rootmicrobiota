

# imports
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, false_discovery_control, wilcoxon
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter

# disable warnings
import warnings
warnings.filterwarnings("ignore")
plt.rcParams.update({
    "font.family": "Helvetica"
})


# load by custom scripts
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
marker_dict = {"native": "s", "non-native": "D"}



# padding of small values
min_val = asv_melted.filter(pl.col("treatment") != "input").filter(pl.col("RA")!= 0)["RA"].min()
max_val = asv_melted.filter(pl.col("treatment") != "input").filter(pl.col("RA")!= 0)["RA"].max()
padding_val = (min_val/max_val)*min_val
asv_melted = asv_melted.with_columns(pl.col("RA").replace(0, padding_val))



# ~~~~~~~~~~~~~~~~~~~ PLOT WITH ALL EVOLVED CANDIDATES ~~~~~~~~~~~~~~~~~~~ #
plot_ls = []
order_ls = []
stats_ls = []
pval_d = {}
pct_increase_d = {}
sizes = []
for syncom in ["at", "lj"]:
    for plant in ["col0", "gifu"]:

        # only perform non-natives
        if ((syncom == "at") and (plant == "col0")) or ((syncom == "lj") and (plant == "gifu")):
            continue

        # get plant strings
        plant_str = "Col-0" if plant == "col0" else "Gifu"
        condition = "native" if ((syncom == "at") and (plant == "col0")) or ((syncom == "lj") and (plant == "gifu")) else "non-native"

        # get the evolved strains subset
        evolved_subset = evolved.filter(pl.col("sc_origin")==syncom).filter(pl.col("plant")==plant)
        strain_ls = evolved_subset["strain"].unique().to_list()
        strain_ls = [strain for strain in strain_ls if strain in sc_strains[syncom].to_list()]


        pval_ls = []
        host_pref_ls = []
        pct_increase_ls = []
        for strain in strain_ls:
            # get median to renormalize data
            ancestral = (asv_melted
                                .filter(pl.col("#Strain")==strain)
                                .filter(pl.col("plant") == plant)
                                .filter(pl.col("treatment") == "at0_vs_lj0")["RA"].median()
                                )
            # normalize
            stats_subest = (asv_melted
                                .filter(pl.col("#Strain")==strain)
                                .filter(pl.col("plant") == plant)
                                .with_columns(pl.col("RA")/ancestral)
                                )

            # perform stats
            pvalue = mannwhitneyu(x=stats_subest.filter(pl.col("treatment") == "at0_vs_lj0")["RA"],
                                  y=stats_subest.filter(pl.col("treatment").str.contains_any([f"{syncom}({plant})"]))["RA"],
                                  alternative="less"
                                ).pvalue
            
            pvalue = wilcoxon(stats_subest.filter(pl.col("treatment").str.contains_any([f"{syncom}({plant})"]))["RA"] - 1, alternative="greater").pvalue
            
            # calculate percent change
            pct_increase = round(((stats_subest.filter(pl.col("treatment").str.contains_any([f"{syncom}({plant})"]))["RA"].median() - stats_subest.filter(pl.col("treatment")=="at0_vs_lj0")["RA"].median()) / stats_subest.filter(pl.col("treatment")=="at0_vs_lj0")["RA"].median()) * 100, 2)

            # save results
            pval_ls.append(pvalue)
            host_pref_ls.append(stats_subest.filter( pl.col("treatment").str.contains_any([f"{syncom}({plant})"]) | (pl.col("treatment") == "at0_vs_lj0" )))
            pct_increase_ls.append(pct_increase)

        # adjust pvalues
        padj_dict = {syncom + "-" + key + "-" + plant: val for key, val in zip(strain_ls, false_discovery_control(pval_ls))}
        pct_increase = {syncom + "-" + key + "-" + plant: val for key, val in zip(strain_ls, pct_increase_ls)}

        # make a stats df to save
        stats_res = (pl.DataFrame({"syncom": syncom, "plant": plant, "strain": strain_ls, "padj": padj_dict.values(), "fitness_c": pct_increase.values()})
                     .filter(pl.col("padj")<0.05)
                     )
        stats_ls.append(stats_res)
        
        # create dict with only significant hits
        pval_d.update({key: val for key, val in padj_dict.items() if val < 0.05})
        pct_increase_d.update({key: val for key, val in pct_increase.items() if padj_dict[key] < 0.05})

        # define order of plotting
        host_pref_df = (pl.concat(host_pref_ls)
                            .with_columns((syncom + "-" + pl.col("#Strain") + "-" + plant).alias("strain"))
                            .filter(pl.col("strain").is_in(pct_increase_d.keys()))
                            .with_columns(pl.lit(condition).alias("condition"))
                       )
        order = host_pref_df.group_by("strain").agg(pl.col("RA").median()).sort("RA", descending=True)["strain"].to_list()
        order_ls += order

        # aggregate data for plotting
        plot_df = (host_pref_df
                    .with_columns(pl.when(pl.col("treatment") == "at0_vs_lj0").then(pl.lit("ancestral")).otherwise(pl.lit("evolved")).alias("strain_type"))
                    .sort(pl.col("strain").cast(pl.Enum(order)), pl.col("strain_type"))
                  )
        plot_ls.append(plot_df)
        sizes.append(len(order))


# ~~~~~~~~~~~~~~~~ MAP MUTATIONS ~~~~~~~~~~~~~~~~ #
# extract mutations for these significant strains from pacbio
pacbio = pl.read_csv("raw_data/pacbio_variants_deepvariant.tsv", separator="\t", schema_overrides={"line":str})
stats_df = (pl.concat(stats_ls)
                .join(pacbio, how="left", on=["syncom", "plant", "strain"])
                .select("strain", "fitness_c", "padj", "chrom", "feature", "pos", "distance",
                        "syncom", "plant", 
                        "geneID", "ko", "alias", "desc",
                        "ref", "alt", "mut_len")
                .rename({"fitness_c": "fitness_change", "pos": "location"})
                .unique()
                .sort("fitness_change", "location", descending=True))
stats_df.write_csv("results/S_adaptive_mutations_pacbio.tsv", separator="\t")
# ~~~~~~~~~~~~~~~~ MAP MUTATIONS ~~~~~~~~~~~~~~~~ #



# aggregate
plot_df = (pl.concat(plot_ls)
            .with_columns(pl.col("plant").alias("evolving_host"))
            .with_columns(pl.when(pl.col("treatment") == "at0_vs_lj0").then(pl.lit("ancestral")).otherwise("plant").alias("plant"))
            .with_columns(pl.col("strain").str.split("-").list[0].alias("syncom"))
            .with_columns((pl.col("syncom") + "-" + pl.col("#Strain") + "-" + pl.col("plant")).alias("strain_x"))
          )

print(plot_df.group_by("evolving_host", "condition", "strain_x").len().sort("len"))
fig, axes = plt.subplots(figsize=(7,4))
p = sns.boxplot(data=plot_df.to_pandas(),
                x="strain_x",
                y="RA",
                hue="evolving_host",
                ax=axes,
                legend=False,
                fill=False,
                showcaps=False,
                palette=color_dict,
                width=0.4,
                **{"showfliers": False}
                )

for c in plot_df["condition"].unique(maintain_order=True):
    sns.stripplot(data=plot_df.filter(pl.col("condition")==c).to_pandas(),
                    x="strain_x",
                    y="RA",
                    hue="evolving_host", 
                    alpha=0.5,
                    marker=marker_dict[c],
                    palette=color_dict,
                    legend=False,
                    ax=axes)

# names on y-axis
axes.set_xticklabels([strain.split("-")[1] + "$^{" + ("At" if strain.split("-")[2] == "col0" else "Lj" if strain.split("-")[2] == "gifu" else "*" )+ "}$" for strain in plot_df["strain_x"].unique(maintain_order=True)], rotation=30, ha="right")
axes.spines[["top", "right"]].set_visible(False)
axes.set_ylim((0, plot_df["RA"].max()*1.35))
axes.yaxis.set_major_formatter(PercentFormatter(xmax=1))

# now add the lines for the comparisons
plot_idx = 0
for strain in plot_df["strain"].unique(maintain_order=True):

    # add the lines
    strain_subset = plot_df.filter(pl.col("strain")==strain)
    strain_max = strain_subset["RA"].max()
    if plot_idx%2 == 0:
        axes.plot((plot_idx,plot_idx+1), 
                  (strain_max*1.15, strain_max*1.15),
                  color="k")
    for strain_x in strain_subset["strain_x"].unique(maintain_order=True):
        substrain_max = strain_subset.filter(pl.col("strain_x") == strain_x)["RA"].max()*1.1
        axes.plot((plot_idx,plot_idx), 
                  (strain_max*1.15, substrain_max),
                  color="k")
        plot_idx+=1
    

    #plot pvalue and pct increase
    print(f"pvalue of {strain}: {pval_d[strain]} , {pretty_exponent(pval_d[strain])}")
    
    axes.text(x=plot_idx-1.5,
              y=strain_max*1.2,
              s=f"{convert_pvalue_to_asterisks(pval_d[strain])}\n+{pct_increase_d[strain]}%",
              ha="center",
              va="bottom")


# set title
axes.set_title("Fitness increase of non-natively evolved isolates")
axes.set_xlabel("")
axes.set_ylabel("Fitness on non-native host")


# custom legend
legend_color = [Patch(label="Evolved host", alpha=0)] + [Patch(facecolor=color, label=get_plant(host), alpha=1) for host, color in color_dict.items()]
lgd = axes.legend(handles=legend_color, 
                  frameon=False,
                  )

for item, label in zip(lgd.legend_handles, lgd.texts):
    if label._text  in ["Evolved host"]:
        width=item.get_window_extent(fig.canvas.get_renderer()).width
        label.set_position((-1.5*width,0))
        label.set_fontsize(12)


# general aesthetics
axes.margins(x=0.025)
fig.tight_layout()
fig.savefig(f"figures/amplicon2_fitness_strain.pdf")
fig.savefig(f"figures/amplicon2_fitness_strain.png", dpi=600)


