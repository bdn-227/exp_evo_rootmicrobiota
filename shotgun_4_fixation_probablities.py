
# ~~~~~~~~~~~~ IMPORTS ~~~~~~~~~~~~ #
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, false_discovery_control, ttest_1samp
import seaborn.objects as so
from matplotlib.patches import Patch
import warnings
import seaborn as sns

# configs
warnings.filterwarnings("ignore")
pl.Config.set_fmt_str_lengths(200)
pl.Config(tbl_rows=100)
plt.rcParams.update({
    "font.family": "Helvetica"
})


# custom functions
from shotgun_analysis_functions import *



# ~~~~~~~~~~~~ LOAD DATA ~~~~~~~~~~~~ #
mapping = load_mapping()
gff = loadGff()


# load the variants
variants_long = (pl.read_csv("raw_data/variants.extended.long.filtered.tsv", 
                             separator="\t", 
                             schema_overrides={"line":str})
                    .filter(pl.col("AF")>0)
                    .filter(pl.col("type")!="upstream")
                    )

# first, reformat variant types
mutation_type_dict = {
        "INTRAGENIC":                           "Intragenic",
        "CODON_CHANGE_PLUS_CODON_INSERTION":    "Insertion",
        "NON_SYNONYMOUS_CODING":                "Non-synonymous",
        "CODON_INSERTION":                      "Insertion",
        "CODON_CHANGE_PLUS_CODON_DELETION":     "Deletion",
        "STOP_GAINED":                          "Nonsense",
        "TRANSCRIPT":                           "Intragenic",
        "FRAME_SHIFT":                          "Frameshift",
        "CODON_DELETION":                       "Deletion",
        'SYNONYMOUS_CODING':                    "Synonymous",
        "START_LOST":                           "Start-lost",
        "INTERGENIC":                           "Intergenic"
    }



# introducing others
mutation_type_dict = {
        "INTRAGENIC":                           "Others",
        "CODON_CHANGE_PLUS_CODON_INSERTION":    "Others",
        "NON_SYNONYMOUS_CODING":                "Non-synonymous",
        "CODON_INSERTION":                      "Others",
        "CODON_CHANGE_PLUS_CODON_DELETION":     "Others",
        "STOP_GAINED":                          "Nonsense",
        "TRANSCRIPT":                           "Others",
        "FRAME_SHIFT":                          "Frameshift",
        "CODON_DELETION":                       "Others",
        'SYNONYMOUS_CODING':                    "Synonymous",
        "START_LOST":                           "Start-lost",
        "INTERGENIC":                           "Intergenic"
    }


variants_long = (variants_long
                 
                        # map the mutation types
                        .with_columns(variants_long["type"].replace(mutation_type_dict).alias("mutation_type"))
                        
                        # add additional column containing information about coding
                        .with_columns(
                                        (pl.when(pl.col("mutation_type") == "Intergenic").then(pl.lit("Non-coding"))
                                           .when(pl.col("mutation_type").is_in(['Frameshift', 'Non-synonymous', 'Start-lost', 'Nonsense', 'Synonymous', "Others"])).then(pl.lit("Coding"))
                                        ).alias("CDS-effect")
                                     )
                        # sort that shit
                        .sort("CDS-effect", "mutation_type")
                        )






# replication of figure 5D from paper
# on level of condition- & strain-level

# load color table
colors = color_dict_strain()
colors["all"] = "k"
alphas = {"native": 0.4, "non-native": 0.8}
colors["all"] = "k"
strain="all"



# another plot
# calculate the fixation probabilities for each condition
grouping = ["treatment", "generation", "line"]
detected_subset = (variants_long
                        .filter(pl.col("AF") > 0)
                        .group_by(grouping)
                        .len("detected")
                        )

fixed_subset = (variants_long
                        .filter(pl.col("AF") >= 0.95)
                        .group_by(grouping)
                        .len("fixed")
                        )

p_fixed = (detected_subset.join(fixed_subset, how="left", on=grouping)
            .fill_null(0)
            .with_columns((pl.col("fixed") / pl.col("detected")).alias("p_fixed"))
            .sort([col for col in ["syncom", "treatment"] if col in grouping])
          )


# test for differences
nat = p_fixed.filter(pl.col("treatment") == "native")["p_fixed"].to_numpy()
nnat = p_fixed.filter(pl.col("treatment") == "non-native")["p_fixed"].to_numpy()
sig = convert_pvalue_to_asterisks(mannwhitneyu(nat, nnat).pvalue)

# figure aesthetic
colors_box = {"native": (0.5, 0.5, 0.5, 1), "non-native": (0, 0, 0, 1)}
colors_points = {"native": (0.5, 0.5, 0.5, 0.3), "non-native": (0, 0, 0, 0.3)}


# make plot
fig, axes = plt.subplots(figsize=(4,3))
print(f"n agg: {p_fixed.shape[0]}")

for e in ["native", "non-native"]:
    sns.stripplot(data=p_fixed.filter(pl.col("treatment")==e).to_pandas(),
                        x="treatment",
                        y="p_fixed",
                        color = colors_points[e],
                        legend=False,
                        ax=axes)
    p = sns.boxplot(data=p_fixed.filter(pl.col("treatment")==e).to_pandas(),
                    x="treatment",
                    y="p_fixed",
                    ax=axes,
                    legend=False,
                    fill=False,
                    color = colors_box[e],
                    width=0.4,
                    showcaps=False,
                    **{"showfliers": False}
                    )



# add statistic results
axes.text(x=.5,
          y=p_fixed["p_fixed"].mean()+0.1,
          s=sig,
          ha="center",
        )

# adjust axis
axes.set_xlabel("")
axes.set_ylabel("$P$[fixed|detected]")
axes.set_ylim((0, 0.6))

# adjust xticks
axes.set_xticklabels([f"{sc.get_text().title()}" for sc in axes.get_xticklabels()], ha="center")

# aesthetics
axes.spines[["top", "right"]].set_visible(False)
axes.set_title(f"Fixation probabilties")

# manage legend
[legend.set_visible(False) for legend in fig.legends]

# just generate the legend
legend_color = [Patch(label="Condition", alpha=0)] + [Patch(facecolor=colors[strain], label=label.title(), alpha=alpha) for label, alpha in alphas.items()]
lgd = axes.legend(handles=legend_color,
                  frameon=False,
                  )
for item, label in zip(lgd.legend_handles, lgd.texts):
    if label._text  in ["Plant", "SynCom", "Generation", "Condition", "Enriched host"]:
        width=item.get_window_extent(fig.canvas.get_renderer()).width
        label.set_position((-1.5*width,0))

# save figures
fig.tight_layout(pad=0)
fig.savefig(f"figures/shotgun_p_aggregated.pdf")
fig.savefig(f"figures/shotgun_p_aggregated.png", dpi=600)













# now strain specific
hyper_mut_idx = ((pl.col("line") == "18") & (pl.col("strain") == "LjRoot33"))
colors = color_dict_strain()
colors_mut = {f"{strain} (mutator)": color for strain, color in colors.items()}
colors_all = {}
colors_all.update(colors)
colors_all.update(colors_mut)



# make the hyper mutator as separate strain
grouping = ["line", "strain"]
detected_subset = (variants_long
                        .with_columns(pl.when(hyper_mut_idx).then(pl.col("strain") + " (mutator)").otherwise("strain").alias("strain"))
                        .filter(pl.col("AF") > 0)
                        .group_by(grouping)
                        .len("detected")
                        )

fixed_subset = (variants_long
                        .with_columns(pl.when(hyper_mut_idx).then(pl.col("strain") + " (mutator)").otherwise("strain").alias("strain"))
                        .filter(pl.col("AF") >= 0.95)
                        .group_by(grouping)
                        .len("fixed")
                        )

p_fixed = (detected_subset.join(fixed_subset, how="left", on=grouping)
            .fill_null(0)
            .with_columns((pl.col("fixed") / pl.col("detected")).alias("p_fixed"))
            #.filter( (pl.col("fixed")!=0) & (pl.col("p_fixed")!=0) ) # remove zero allele frequencies
          )


# test for difference in hypermutator
mutator_ls = list(set([strain for strain in p_fixed["strain"].to_list() if "mutator" in strain]))
pval_ls = []
strain_ls = []
for mutator in mutator_ls:
    strain = mutator.split(" ")[0]
    mutator_set = p_fixed.filter(pl.col("strain") == mutator)
    strain_set = p_fixed.filter(pl.col("strain") == strain)
    pval_ls.append(ttest_1samp(a=strain_set["p_fixed"], popmean=mutator_set["p_fixed"].mean()).pvalue)
    strain_ls.append(strain)
sig_vec = np.array([convert_pvalue_to_asterisks(padj) for padj in  false_discovery_control(pval_ls)])
mutator_ls = np.array(mutator_ls)[sig_vec!="NS"]
strain_ls = np.array(strain_ls)[sig_vec!="NS"]
sig_vec = np.array(sig_vec)[sig_vec!="NS"]


#order
order = p_fixed.group_by("strain").agg(pl.col("p_fixed").mean()).sort("p_fixed", descending=True)["strain"]

# make plot
print(f"n_strain: {p_fixed.shape[0]}")
fig, axes = plt.subplots(figsize=(8, 4))
(so.Plot(data=p_fixed.sort(pl.col("strain").cast(pl.Enum(order))).to_pandas(),
         x="strain",
         y="p_fixed",
         color="strain",
        )
        
        # add the artists
        .add(so.Dot(edgewidth=0, pointsize=3, alpha=0.6), so.Jitter(), legend=False)
        .add(so.Bar(edgewidth=0, alpha=0.6), so.Agg(func="mean"), legend=False)
        .add(so.Range(alpha=1), so.Est(errorbar=("ci", 95), func="mean"), legend=False)

        .scale(color=colors_all)

        .on(axes)
        ).plot()


# plot N's
for idx, strain in enumerate(order):
    axes.text(x=idx,
              y=1.05,
              s=f'$n$={p_fixed.filter(pl.col("strain")==strain)["detected"].max()}',
              ha="center",
              fontsize=4,
             )


# adjust axis
axes.set_xlabel("")
axes.set_ylabel("$P$[fixed|detected]")
axes.set_ylim((0, 1.1))

# adjust xticks
names_df = load_colors()
strain_name_dict = {strain: f"{strain}" for strain, family in zip(names_df["strain"], names_df["family"])}
mutator_name_dict = {f"{strain} (mutator)": f"{strain} (mutator)" for strain, family in zip(names_df["strain"], names_df["family"])}
names_all = {}
names_all.update(strain_name_dict)
names_all.update(mutator_name_dict)


xlim = axes.get_xlim()

# aesthetics
axes.spines[["top", "right"]].set_visible(False)
axes.set_title(f"Strain specific fixation probabilties")


# add statistic results for mutator
for strain, mutator, sig in zip(strain_ls.tolist(), mutator_ls.tolist(), sig_vec.tolist()):
    xcoord = [label.get_position()[0] for label in axes.get_xticklabels() if label.get_text() in (strain, mutator)]
    y_pos = p_fixed.filter(pl.col("strain").is_in([strain, mutator]))["p_fixed"].max()
    strain_1 = [label.get_text() for label in axes.get_xticklabels() if label.get_position()[0]==xcoord[0]]
    strain_2 = [label.get_text() for label in axes.get_xticklabels() if label.get_position()[0]==xcoord[1]]
    axes.plot(xcoord, [y_pos+0.075]*2, color="k", linewidth=1)
    axes.plot([xcoord[0], xcoord[0]], [p_fixed.filter(pl.col("strain").is_in(strain_1))["p_fixed"].max()+0.05, y_pos+0.075], color="k", linewidth=1)
    axes.plot([xcoord[1], xcoord[1]], [p_fixed.filter(pl.col("strain").is_in(strain_2))["p_fixed"].max()+0.05, y_pos+0.075], color="k", linewidth=1)
    axes.text(np.mean(xcoord), y=y_pos+0.1, s=sig, ha="center")

axes.set_xticklabels([names_all[label.get_text()] for label in axes.get_xticklabels()], rotation=45, ha="right")

# save figures
axes.set_xlim(xlim)
fig.tight_layout(pad=0.75)
fig.savefig(f"figures/shotgun_p_strain.pdf")
fig.savefig(f"figures/shotgun_p_strain.png", dpi=600)


# get the values
def get_family_probability(family):
    strains = load_colors().filter(pl.col("family")==family)["strain"].unique().to_list()
    return p_fixed.filter(pl.col("strain").is_in(strains))["p_fixed"].max()*100

get_family_probability("Flavobacteriaceae")


# compare ljRoot33 vs the hyper mutator
print(f'fixation probabilites: LjRoot33: {round(p_fixed.filter(pl.col("strain") == "LjRoot33")["p_fixed"].mean()*100, 2)}%')
print(f'fixation probabilites: LjRoot33(mutator) {round(p_fixed.filter(pl.col("strain") == "LjRoot33 (mutator)")["p_fixed"].mean()*100, 2)}%')


