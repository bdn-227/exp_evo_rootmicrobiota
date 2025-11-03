
# ~~~~~~~~~~~~ IMPORTS ~~~~~~~~~~~~ #
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import mannwhitneyu, false_discovery_control
import warnings
warnings.filterwarnings("ignore")

# configs
#mpl.use("TkAgg")
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
                        .with_columns(variants_long["type"].replace(mutation_type_dict).alias("mutation_type"))
                        .with_columns(
                                        (pl.when(pl.col("mutation_type") == "Intergenic").then(pl.lit("Non-coding"))
                                           .when(pl.col("mutation_type").is_in(['Frameshift', 'Non-synonymous', 'Start-lost', 'Nonsense', 'Synonymous', "Others"])).then(pl.lit("Coding"))
                                        ).alias("CDS-effect")
                                     )
                        .sort("CDS-effect", "mutation_type")
                        )


# get the genomic coding/non-coding region
sc_size = (pl.read_csv("mapping/syncom_size.txt", separator="\t")
                .filter(pl.col("syncom") == "mixed")
                .with_columns((pl.col("sum_len_cds") / pl.col("sum_length_genome")).alias("Coding"))
                .with_columns((1-pl.col("Coding")).alias("Non-coding"))
           )


# order of things
order = ['generation', 'Frameshift', 'Non-synonymous', 'Nonsense', 'Start-lost', 'Synonymous', 'Others', 'Intergenic']
order_coding = ["generation", "Coding", "Non-coding"]


# define color map
color_map = {mutation_type: mpl.colormaps["Set3"](index) for index, mutation_type in enumerate(variants_long["mutation_type"].unique(maintain_order=True))}
color_map_coding = {"Non-coding": "grey", "Coding": "black"}




ratios_syn = []
ratios_coding = []
coding_fraction_d = {}
for idx, cutoff in enumerate([0, 0.95]):

    fig, axes = plt.subplots(nrows=1,
                             ncols=1,
                             figsize=(5,3)
                             )

    # second set of data for coding block
    coding_info = (variants_long
                        .filter(pl.col("AF") >= cutoff) 
                        .group_by("generation", "CDS-effect")
                        .len()
                        .pivot(index="generation", on="CDS-effect")
                        .fill_null(0)
                        .sort("generation")  
                        .to_pandas()
                        .set_index("generation") )
    

    coding_fraction_d[cutoff] = coding_info["Coding"].sum()/(coding_info["Coding"].sum() + coding_info["Non-coding"].sum())
    
    # get the last generation and rescale the genomic background
    gens = coding_info.index
    total = coding_info.iloc[-1,:].sum()
    background = (sc_size.select("Coding", "Non-coding") * total).with_columns(pl.lit("Genome").alias("generation")).to_pandas().set_index("generation")

    # make stacked barplot using pandas
    plot_data = (variants_long
                    .filter(pl.col("AF") >= cutoff) 
                    .group_by("generation", "mutation_type")
                    .len()
                    .pivot(index="generation", on="mutation_type")
                    .fill_null(0)
                    .sort("generation")
                    .to_pandas()
                    .set_index("generation") )
    

    # perform statistics on the ratio of coding / non-coding OR synonymous vs all others
    ratio = np.array(plot_data["Non-synonymous"] / plot_data["Synonymous"])
    ratio = ratio[~np.isinf(ratio)]
    ratio = ratio[~np.isnan(ratio)]
    ratios_syn.append(ratio)

    ratio = np.array(coding_info["Coding"][gens] / coding_info["Non-coding"][gens])
    ratio = ratio[~np.isinf(ratio)]
    ratio = ratio[~np.isnan(ratio)]
    ratios_coding.append(ratio)

    # define bottom
    coding_info = coding_info[[col for col in order_coding if col in coding_info.columns]]
    bottom = np.zeros(coding_info.iloc[:,0].shape)

    # stack them myself so I can add edge colors
    for coding_type in coding_info.columns:
        plt_subset = coding_info[coding_type]#[coding_info[coding_type] > 0]
        axes.bar(
                 x=plt_subset.index.astype(str),
                 height=plt_subset,
                 color=color_map_coding[coding_type],
                 edgecolor=color_map_coding[coding_type],
                 linewidth=.9,
                 alpha=1,
                 bottom=bottom[np.arange(plt_subset.index.shape[0])],
                 width=0.4,
                 align="center",
                 )
        bottom[np.arange(plt_subset.index.shape[0])] += plt_subset

    # define bottom
    plot_data = plot_data[[col for col in order if col in plot_data.columns]]
    bottom = np.zeros(plot_data.iloc[:,0].shape)

    # stack them myself so I can add edge colors
    for mutation_type in plot_data.columns:
        plt_subset = plot_data[mutation_type][plot_data[mutation_type] > 0]
        axes.bar(
                 x=plt_subset.index.astype(str),
                 height=plt_subset,
                 color=color_map[mutation_type],
                 edgecolor=color_map[mutation_type],
                 linewidth=1,
                 alpha=1,
                 bottom=bottom[plt_subset.index - 1],
                 width=-0.4,
                 align="edge",
                 )
        bottom[plt_subset.index - 1] += plt_subset

    # set labels
    axes.set_ylabel("Number of mutations")
    axes.set_xlabel("Plant cycle, $t$")

    # x-limits
    axes.margins(x=0.01)
    axes.spines[["top", "right"]].set_visible(False)

    # set xticks
    #axes.xaxis.get_major_ticks()[-2].set_visible(False)
    xticklabels=[]
    for text in axes.get_xticklabels():
        if text.get_text() not in (list(gens.astype(str)) + ["Genome"]):
            text.set_text("")
        xticklabels.append(text)
    axes.set_xticklabels(xticklabels)


    # add titles
    if idx==0:
        fig.suptitle("All mutations")
    else:
        fig.suptitle(f"Fixed mutations (AF ≥ {cutoff})")

    fig.tight_layout(pad=0)
    fig.savefig(f"figures/shotgun_variant_types_{cutoff}.png", dpi=600)
    fig.savefig(f"figures/shotgun_variant_types_{cutoff}.pdf")



# perform statstics
pval_syn = mannwhitneyu(ratios_syn[0], ratios_syn[1]).pvalue
pval_coding = mannwhitneyu(ratios_coding[0], ratios_coding[1]).pvalue
padj_syn, padj_coding = false_discovery_control(np.array([pval_syn, pval_coding]))
title = f"Difference in variant types between all and fixed variants\nNon-coding vs Coding $P$-value={padj_coding}\nNon-synonymous vs Synonymous $P$-value={padj_syn}"
print(title)


def fold_increase(old, new):
    if old == 0:
        raise ValueError("Old value cannot be zero when calculating fold increase.")
    return new / old

# print coding fractions
coding_fraction_d.keys()
print(f"fraction of genes affecting coding regions for all variants: {round(coding_fraction_d[0]*100, 2)}%")
print(f"fraction of genes affecting coding regions for fixed variants: {round(coding_fraction_d[0.95]*100, 2)}%")


# calulcate fold-change
print(f"increase in ratio of non-syn to syn: {round(fold_increase(np.mean(ratios_syn[0]), np.mean(ratios_syn[1])), 2)}")



# save legend separately
fig, axes = plt.subplots(nrows=1,
                         ncols=1, 
                         figsize=(1.8,2.2), 
                        )

axes.xaxis.set_visible(False)
axes.yaxis.set_visible(False)
axes.spines[["top", "bottom", "left", "right"]].set_visible(False)

# add legend
class Handler(object):
    def __init__(self, color1, color2, only_second=False):
        self.color1=color1
        self.color2=color2
        self.only_second=only_second
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        if self.only_second:
            patch = plt.Rectangle([x0+width, y0], width/4., height, facecolor=self.color2, edgecolor='k', transform=handlebox.get_transform())
            handlebox.add_artist(patch)
        else:
            patch = plt.Rectangle([x0, y0], width, height, facecolor=self.color1, edgecolor='k', transform=handlebox.get_transform())
            patch2 = plt.Rectangle([x0+width, y0], width/4., height, facecolor=self.color2, edgecolor='k', transform=handlebox.get_transform())
            handlebox.add_artist(patch)
            handlebox.add_artist(patch2)
        return patch


# add cds mappings to legend
color_map.update(color_map_coding)
order_vec = np.array(order[1:])[::-1]
order_vec = list(order_vec) + list(color_map_coding.keys())

handles = [plt.Rectangle((0,0),1,1) for _ in order_vec]
hmap = dict(zip(handles, [Handler(color_map[label], "grey" if label in ["Intergenic","Non-coding"] else "black", True if label in ["Coding","Non-coding"] else False) for label in order_vec]))

lgd = axes.legend(handles=handles, 
            labels=order_vec, 
            handler_map=hmap, 
            loc='center', 
            bbox_to_anchor=(0.5, 0.5), 
            title="Type of Mutation",
            frameon=False)

for item, label in zip(lgd.legend_handles, lgd.texts):
    if label._text  in ["Type of Mutation"]:
        width=item.get_window_extent(fig.canvas.get_renderer()).width
        label.set_position((-2.0*width,0))


fig.tight_layout()
fig.savefig(f"figures/shotgun_variant_types_legend.png", dpi=600, bbox_inches="tight")
fig.savefig(f"figures/shotgun_variant_types_legend.pdf")


