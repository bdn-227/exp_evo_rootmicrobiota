
# import modules
import polars as pl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings("ignore")


# ~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~ #
from shotgun_analysis_functions import *
plt.rcParams.update({
    "font.family": "Helvetica"
})


# ~~~~~~~~~~~~ MAPPING ~~~~~~~~~~~~ #
mapping = load_mapping()


# and import the inputs
file="raw_data/variants.extended.long.filtered.tsv"
variants_long = (pl.read_csv(file, 
                             separator="\t",
                             schema_overrides={"line": str}
                            )
                            .filter(pl.col("type")!="upstream")
                            .filter(pl.col("af_type")=="original")
                            .filter(pl.col("line").is_in(["1","2","3","4","5",
                                                          "6.2","7.2","8.2","9","10",
                                                          "11","12","13","14","15",
                                                          "16","17","18","19","20"]))
                            .select("ID", "plant", "line", "syncom", "generation", "AF")
                            )


# load ancestrals
file=f"raw_data/variants.gen_0.filtered.long.tsv"
ancestral_long = (pl.read_csv(file, separator="\t")
                            .filter(pl.col("type")!="upstream")
                            .select(variants_long.columns)
                            )


# plot based on allele frequence values
variants = (pl.concat([variants_long, ancestral_long], how="vertical")
                    .with_columns( (pl.col("line") + "_g" + pl.col("generation").cast(pl.String)).alias("line_gen") ) 
                    .select("AF", "line_gen", "ID", "plant", "syncom", "generation") 
                    .pivot(on="ID", values="AF")
                    .with_columns( pl.when(pl.col("plant") == "col0").then(1).when(pl.col("plant") == "gifu").then(2).otherwise(3).alias("plant_encoded") )
                    .with_columns( pl.when(pl.col("syncom") == "at").then(1).when( pl.col("syncom") =="lj").then(2).otherwise(3).alias("syncom_encoded") )
                    .fill_null(0)
                    )



points_ls = []
variance = []
for syncom in ["at", "lj"]:

    # split the dataframe
    variants_wide = variants.filter(pl.col("syncom") == syncom).to_pandas().set_index("line_gen")
    variants_matrix = variants_wide.drop(columns=["plant", "syncom", "generation", "plant_encoded", "syncom_encoded"])
    variants_metadata = variants_wide[["generation", "plant_encoded", "syncom_encoded", "plant", "syncom"]]

    # calculate distance matrix
    full_dm = squareform(pdist(variants_matrix.to_numpy(), metric="euclidean"))
    full_dm = DistanceMatrix(full_dm, ids=variants_matrix.index)
    
    # perform pcoa
    pcoa_res = pcoa(full_dm.to_data_frame())
    
    # calcluate variance of axis and total
    variance_explained = pcoa_res.proportion_explained*100
    total_var = pcoa_res.eigvals.sum()/(pcoa_res.eigvals.sum() + pcoa_res.eigvals.sum())
    total_var = f"{round(total_var*100, 2)}%"
    variance.append(total_var)
    

    points = pl.DataFrame(pcoa_res.samples)
    points.columns = [f"PC{i+1}" for i in range(len(points.columns))]
    points = pl.concat([points, pl.DataFrame(variants_metadata.reset_index())], how="horizontal").with_columns(pl.lit(variance_explained[0]).alias("var_x")).with_columns(pl.lit(variance_explained[1]).alias("var_y"))
    points_ls.append(points)


# concat stuff
points_df = (pl.concat(points_ls, how="diagonal_relaxed").with_columns( (pl.when( (pl.col("plant") == "col0") & (pl.col("syncom") == "at") ).then(pl.lit("native"))
                                                                           .when( (pl.col("plant") == "gifu") & (pl.col("syncom") == "lj") ).then(pl.lit("native"))
                                                                           .when( pl.col("line_gen").str.contains("zero") ).then(pl.lit("zero"))
                                                                           .when( pl.col("plant") == "none" ).then(pl.lit("ancestral"))
                                                                           .otherwise(pl.lit("non-native"))).alias("treatment")
                                                                        ) )

color_dict = {"col0": "red", "gifu": "blue", "none": "black"}
markers = {"ancestral": "o", "native": "s", "non-native": "D"} #"zero": "*",

size_dict = {}
size_dict["Zero"] = 50
size_dict["0"] = 50
for generation in range(1,17):
    size_dict[str(generation)] = (20+(generation/1.5)**2)





for idx, sc in enumerate(["at", "lj"]):

    # make figure
    fig, axes = plt.subplots(figsize=(5,5))
    plt_subset = points_df.filter(pl.col("syncom") == sc).with_columns( pl.col("generation").cast(str) ).to_pandas()
    print(f"total mutations found in {sc}: " + str(plt_subset.shape[0]))
    sns.scatterplot(data=plt_subset,
                    x="PC1",
                    y="PC2",
                    alpha=0.6,
                    
                    hue="plant",
                    palette=color_dict,

                    markers=markers,
                    style="treatment",
                    size="generation",
                    sizes=size_dict,
                    linewidth=0,


                    legend=False,
                    ax=axes)
    
    axes.set_xlabel(f"PC1: {np.round(plt_subset['var_x'][0], 2)}%")
    axes.set_ylabel(f"PC2: {np.round(plt_subset['var_y'][0], 2)}%")
    axes.set_title(f"Variant profiles for ${(sc.title())}$-SC")
    axes.spines[["top", "right"]].set_visible(False)
    

    # save the figure
    fig.tight_layout()
    fig.savefig(f"figures/shotgun_variant_pcoa_AF_{sc}.pdf")
    fig.savefig(f"figures/shotgun_variant_pcoa_AF_{sc}.png", dpi=600)




# make custom legend
size_dict2 = {}
size_dict2["Zero/Ancestral"] = size_dict["Zero"]
size_dict.pop("Zero")
size_dict.pop("0")
size_dict2.update(size_dict)
legend_color = [Patch(label="Plant", alpha=0)] + [Patch(facecolor=color, label=get_plant_label(label), alpha=0.6) for label, color in color_dict.items()]
legend_shape = [Patch(label="Condition", alpha=0)] + [Line2D([0], [0], marker=marker, markersize=8 if label !="zero" else 14, color='w', label=label.title(), markerfacecolor='black') for label, marker in markers.items()]
legend_size = [Patch(label="Plant cycle", alpha=0)] + [Line2D([0], [0], marker="o", color='w', label=label , markerfacecolor='black', markersize=np.sqrt(size) ) for label, size in size_dict2.items() if label not in np.arange(1,16,2).astype(str)]



# just generate the legend
fig, axes = plt.subplots(figsize=(2,4.1), 
                        )

lgd = axes.legend(handles=legend_color + legend_shape + legend_size, 
                  loc='center', 
                  bbox_to_anchor=(0.5, 0.5),
                  frameon=False,
                  )
axes.xaxis.set_visible(False)
axes.yaxis.set_visible(False)
axes.spines[["top", "bottom", "left", "right"]].set_visible(False)

for item, label in zip(lgd.legend_handles, lgd.texts):
    if label._text  in ["Plant", "SynCom", "Plant cycle", "Condition"]:
        width=item.get_window_extent(fig.canvas.get_renderer()).width
        label.set_position((-1.5*width,0))
        label.set_fontsize(12)

fig.tight_layout()
fig.savefig(f"figures/shotgun_variant_pcoa_AF_legend.png", dpi=600)
fig.savefig(f"figures/shotgun_variant_pcoa_AF_legend.pdf")











